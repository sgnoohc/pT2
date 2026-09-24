import uproot
import awkward as ak
import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import IterableDataset, DataLoader
from sklearn.metrics import roc_auc_score, roc_curve, precision_recall_curve
from sklearn.calibration import calibration_curve
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import LogNorm
import argparse
import glob
import os
import random
import json
import torch.onnx
from tqdm import tqdm
from collections import deque

# ---------------------------
# Args
# ---------------------------
parser = argparse.ArgumentParser(
    description="Train the pT2 real-vs-fake NN on pt2_training_data.root files from `pt2 process -r`")
parser.add_argument("--data",        type=str, nargs="+", default=["data_new/"],
                    help="Training files, directories of .root files, or glob patterns")
parser.add_argument("--output_dir",  type=str, default="outputs/")
parser.add_argument("--val_frac",    type=float, default=0.1, help="Fraction of events used for validation")
parser.add_argument("--test_frac",   type=float, default=0.1, help="Fraction of events used for testing")
parser.add_argument("--epochs",      type=int, default=30)
parser.add_argument("--fast_dev_run", action="store_true",
                    help="Quick check: 2 epochs, at most 20 batches per epoch")
parser.add_argument("--skip_shap",   action="store_true",
                    help="Skip SHAP (slow on large test sets; needs the shap package)")
args = parser.parse_args()

os.makedirs(args.output_dir, exist_ok=True)
plot_dir = os.path.join(args.output_dir, "plots")
os.makedirs(plot_dir, exist_ok=True)

DEVICE     = "cuda" if torch.cuda.is_available() else "cpu"
BATCH_SIZE = 32768
EPOCHS     = 2 if args.fast_dev_run else args.epochs
MAX_BATCHES = 20 if args.fast_dev_run else None   # per pass over a split

# Rows read from a file at a time (keeps memory bounded for large files)
READ_STEP = 2_000_000

# Rows are split into train / val / test by event, so all pT2s of one event land in the same split
SPLIT_SEED = 12345

# True global ratio: ~1.2B fakes / ~3M reals
POS_WEIGHT = 400.0

# Score gap penalty: weight of auxiliary loss pushing real scores above this threshold
GAP_PENALTY_WEIGHT     = 0.5   # relative weight vs BCE loss
GAP_PENALTY_THRESHOLD  = 0.9   # penalise real tracks scoring below this

# Log-transform epsilon
LOG_EPS = 1e-6

# ---------------------------
# Features
# v7 change: added log_abs_md0_dxy (single highest-SHAP feature, log scale)
# ---------------------------
FEATURE_NAMES = [
    "ls_pt", "ls_eta", "ls_sin_phi", "ls_cos_phi",
    "pls_pt", "pls_eta", "pls_sin_phi", "pls_cos_phi",
    "pls_charge", "pls_nhit",
    "pt2_delta_pt", "pt2_delta_eta", "pt2_delta_phi", "pt2_delta_R",
    "pt2_md0_dxy", "pt2_md0_dz",
    "pt2_md1_dxy", "pt2_md1_dz",
    "pt2_md0_rz",  "pt2_md1_rz",
    "log_abs_md0_dxy",   # NEW
]

# ---------------------------
# Plot style
# ---------------------------
STYLE = {
    "real_color":   "#1D9E75",
    "fake_color":   "#D85A30",
    "accent":       "#378ADD",
    "neutral":      "#888780",
    "threshold_99": "#E24B4A",
    "threshold_95": "#BA7517",
    "threshold_90": "#7F77DD",
    "threshold_999":"#0F6E56",
    "bg":           "#F8F8F6",
    "grid_alpha":   0.25,
}

def style_ax(ax, title=None, xlabel=None, ylabel=None):
    ax.set_facecolor(STYLE["bg"])
    ax.grid(True, alpha=STYLE["grid_alpha"], linewidth=0.5)
    ax.spines[["top","right"]].set_visible(False)
    if title:  ax.set_title(title,  fontsize=11, fontweight="normal", pad=6)
    if xlabel: ax.set_xlabel(xlabel, fontsize=10)
    if ylabel: ax.set_ylabel(ylabel, fontsize=10)
    ax.tick_params(labelsize=9)

def savefig(name):
    plt.savefig(os.path.join(plot_dir, name), dpi=150, bbox_inches="tight",
                facecolor="white")
    plt.close()

# ---------------------------
# Dataset
# ---------------------------
BRANCHES = [
    "event_idx", "is_real",
    "ls_pt", "ls_eta", "ls_phi", "pls_pt", "pls_eta", "pls_phi", "pls_charge", "pls_nhit",
    "pt2_delta_pt", "pt2_delta_eta", "pt2_delta_phi", "pt2_delta_R",
    "pt2_md0_dxy", "pt2_md0_dz", "pt2_md1_dxy", "pt2_md1_dz", "pt2_md0_rz", "pt2_md1_rz",
]

def event_split(event_idx):
    """0 = train, 1 = val, 2 = test; deterministic per event index."""
    u = (event_idx.astype(np.uint64) * np.uint64(2654435761) + np.uint64(SPLIT_SEED)) % np.uint64(1_000_003)
    u = u.astype(np.float64) / 1_000_003
    return np.where(u < args.test_frac, 2, np.where(u < args.test_frac + args.val_frac, 1, 0))

def make_features(d):
    """Same 21 inputs, in the same order, as Pt2Scorer::features() in src/pt2_scorer.cc."""
    ls_phi, pls_phi = d["ls_phi"], d["pls_phi"]
    md0_dxy = d["pt2_md0_dxy"]
    x = np.stack([
        d["ls_pt"],  d["ls_eta"],
        np.sin(ls_phi), np.cos(ls_phi),
        d["pls_pt"], d["pls_eta"],
        np.sin(pls_phi), np.cos(pls_phi),
        d["pls_charge"], d["pls_nhit"],
        d["pt2_delta_pt"],  d["pt2_delta_eta"],
        d["pt2_delta_phi"], d["pt2_delta_R"],
        md0_dxy,            d["pt2_md0_dz"],
        d["pt2_md1_dxy"],   d["pt2_md1_dz"],
        d["pt2_md0_rz"],    d["pt2_md1_rz"],
        np.log(np.abs(md0_dxy) + np.float32(LOG_EPS)),
    ], axis=1).astype(np.float32)
    # Pt2Scorer replaces non-finite inputs by 0 before normalizing; do the same
    return np.nan_to_num(x, nan=0.0, posinf=0.0, neginf=0.0)

class Pt2Dataset(IterableDataset):
    """Yields (features, is_real) batches of one split, shuffled within each read chunk."""
    def __init__(self, files, split):
        self.files = files
        self.split = split

    def __iter__(self):
        n_batches = 0
        for d in uproot.iterate([f + ":tree" for f in self.files], BRANCHES,
                                step_size=READ_STEP, library="np"):
            keep = event_split(d["event_idx"]) == self.split
            if not keep.any():
                continue
            is_real = d["is_real"][keep].astype(np.float32)
            features = make_features({k: v[keep] for k, v in d.items()})

            idx = np.random.permutation(len(is_real))
            features, is_real = features[idx], is_real[idx]

            for i in range(0, len(is_real), BATCH_SIZE):
                yield (torch.from_numpy(features[i:i+BATCH_SIZE]),
                       torch.from_numpy(is_real[i:i+BATCH_SIZE]))
                n_batches += 1
                if MAX_BATCHES is not None and n_batches >= MAX_BATCHES:
                    return

# ---------------------------
# Model — with residual connection
# ---------------------------
class Model(nn.Module):
    """
    v7 architecture:
      - Block 1: Linear(input→256) + BN + ReLU + Dropout(0.2)
      - Block 2: Linear(256→128)   + BN + ReLU + Dropout(0.1)
      - Residual projection: Linear(input→128, no bias) added after Block 2
      - Block 3: Linear(128→64)    + ReLU
      - Head:    Linear(64→1)

    The residual projection gives gradients a direct path back to the input,
    which helps sharpen score separation at high efficiency working points.
    """
    def __init__(self, input_dim):
        super().__init__()

        self.block1 = nn.Sequential(
            nn.Linear(input_dim, 256),
            nn.BatchNorm1d(256),
            nn.ReLU(),
            nn.Dropout(0.2),        # reduced from 0.3
        )
        self.block2 = nn.Sequential(
            nn.Linear(256, 128),
            nn.BatchNorm1d(128),
            nn.ReLU(),
            nn.Dropout(0.1),        # reduced from 0.2
        )
        # Residual projection: maps raw input directly to 128-dim space
        self.residual_proj = nn.Linear(input_dim, 128, bias=False)

        self.block3 = nn.Sequential(
            nn.Linear(128, 64),
            nn.ReLU(),
        )
        self.head = nn.Linear(64, 1)

    def forward(self, x):
        out = self.block1(x)
        out = self.block2(out) + self.residual_proj(x)   # residual add
        out = self.block3(out)
        return self.head(out).squeeze(-1)


# ---------------------------
# SHAP shim
# ---------------------------
class ShapShim(nn.Module):
    def __init__(self, model):
        super().__init__()
        self.model = model

    def forward(self, x):
        return self.model(x).unsqueeze(-1)


# ---------------------------
# Loss: BCE + score gap penalty
# ---------------------------
def gap_penalty(logits, y):
    """
    Auxiliary loss that penalises real tracks (y==1) whose sigmoid score
    falls below GAP_PENALTY_THRESHOLD.  This pushes all real tracks toward
    score=1, closing the gap in the score distribution that causes the
    99% threshold to fall into noisy fake territory.
    """
    scores     = torch.sigmoid(logits)
    real_mask  = y == 1
    if real_mask.sum() == 0:
        return torch.tensor(0.0, device=logits.device)
    real_scores = scores[real_mask]
    # Hinge-style: only penalise scores below threshold
    penalty = torch.clamp(GAP_PENALTY_THRESHOLD - real_scores, min=0.0)
    return penalty.mean()


def combined_loss(criterion, logits, y):
    bce     = criterion(logits, y)
    penalty = gap_penalty(logits, y)
    return bce + GAP_PENALTY_WEIGHT * penalty


# ---------------------------
# Train / Eval
# ---------------------------
def train_epoch(model, loader, optimizer, criterion):
    model.train()
    losses = []
    for x, y in loader:
        x, y = x.to(DEVICE), y.to(DEVICE)
        optimizer.zero_grad()
        loss = combined_loss(criterion, model(x), y)
        loss.backward()
        optimizer.step()
        losses.append(loss.item())
    return np.mean(losses)

def evaluate(model, loader):
    model.eval()
    ys, preds = [], []
    with torch.no_grad():
        for x, y in loader:
            x = x.to(DEVICE)
            out = torch.sigmoid(model(x))
            ys.append(y.cpu().numpy())
            preds.append(out.cpu().numpy())
    ys    = np.concatenate(ys)
    preds = np.concatenate(preds)
    return roc_auc_score(ys, preds), ys, preds

def fake_rate_at_efficiency(y_true, y_pred, target_eff):
    real_scores = y_pred[y_true == 1]
    fake_scores = y_pred[y_true == 0]
    threshold   = np.quantile(real_scores, 1 - target_eff)
    fake_eff    = (fake_scores > threshold).mean()
    return fake_eff, threshold

# ---------------------------
# Load files
# ---------------------------
files = []
for pattern in args.data:
    if os.path.isdir(pattern):
        files += glob.glob(os.path.join(pattern, "*.root"))
    else:
        files += glob.glob(pattern)
files = sorted(set(files))
if not files:
    raise SystemExit(f"No .root files found in {args.data}")

print(f"{len(files)} file(s); event split train/val/test = "
      f"{1 - args.val_frac - args.test_frac:.2f}/{args.val_frac:.2f}/{args.test_frac:.2f}")
train_files = val_files = test_files = files
TRAIN, VAL, TEST = 0, 1, 2

# Every split needs at least one event, and training needs real pT2s
split_counts = np.zeros((3, 2), dtype=np.int64)   # [split][events, reals]
for f in files:
    events = set()
    for d in uproot.iterate(f + ":tree", ["event_idx", "is_real"], step_size=READ_STEP, library="np"):
        events.update(np.unique(d["event_idx"]).tolist())
        np.add.at(split_counts[:, 1], event_split(d["event_idx"]), d["is_real"].astype(np.int64))
    np.add.at(split_counts[:, 0], event_split(np.array(sorted(events), dtype=np.int64)), 1)
for name, k in [("train", TRAIN), ("val", VAL), ("test", TEST)]:
    print(f"  {name:5s}: {split_counts[k, 0]} events, {split_counts[k, 1]} real pT2s")
if (split_counts[:, 0] == 0).any() or split_counts[TRAIN, 1] == 0:
    raise SystemExit("Not enough events for a train/val/test split; "
                     "produce more with `pt2 process -r -n <N>` or adjust --val_frac/--test_frac")

# ---------------------------
# Normalization — dedicated loader, not reused
# ---------------------------
print("Computing normalization...")
norm_loader = DataLoader(Pt2Dataset(train_files, TRAIN), batch_size=None)
all_feats = []
for i, (x, _) in enumerate(norm_loader):
    all_feats.append(x.numpy())
    if i > 50: break
del norm_loader

all_feats = np.concatenate(all_feats)
mean = all_feats.mean(axis=0)
std  = all_feats.std(axis=0) + 1e-6

np.save(os.path.join(args.output_dir, "mean.npy"), mean)
np.save(os.path.join(args.output_dir, "std.npy"),  std)

def normalize(x):
    return (x - mean) / std

class NormWrapper:
    def __init__(self, loader):
        self.loader = loader
    def __iter__(self):
        for x, y in self.loader:
            yield torch.tensor(normalize(x.numpy()), dtype=torch.float32), y

train_loader = NormWrapper(DataLoader(Pt2Dataset(train_files, TRAIN), batch_size=None))
val_loader   = NormWrapper(DataLoader(Pt2Dataset(val_files,   VAL),   batch_size=None))
test_loader  = NormWrapper(DataLoader(Pt2Dataset(test_files,  TEST),  batch_size=None))

# ---------------------------
# Model setup
# ---------------------------
input_dim = len(mean)
assert input_dim == len(FEATURE_NAMES), \
    f"Feature mismatch: mean has {input_dim}, FEATURE_NAMES has {len(FEATURE_NAMES)}"
print(f"Input dim: {input_dim}")

model     = Model(input_dim).to(DEVICE)
optimizer = torch.optim.Adam(model.parameters(), lr=1e-3)
criterion = nn.BCEWithLogitsLoss(
    pos_weight=torch.tensor([POS_WEIGHT]).to(DEVICE)
)
scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(
    optimizer, T_max=EPOCHS, eta_min=1e-5
)

# ---------------------------
# Training loop
# v7: early stopping on val AUC (stable), patience=8
#     fake rate still logged for visibility
# ---------------------------
train_losses, val_aucs, val_fake_rates, lrs = [], [], [], []
best_val_auc   = -float("inf")
patience       = 8
min_delta      = 1e-5
no_improve     = 0
fake_rate_window = deque(maxlen=3)

for epoch in range(EPOCHS):
    loss = train_epoch(model, train_loader, optimizer, criterion)
    auc, y_val, p_val = evaluate(model, val_loader)
    scheduler.step()

    fake_rate, _ = fake_rate_at_efficiency(y_val, p_val, target_eff=0.99)
    current_lr   = scheduler.get_last_lr()[0]

    fake_rate_window.append(fake_rate)
    rolling_fr = float(np.mean(fake_rate_window))

    train_losses.append(loss)
    val_aucs.append(auc)
    val_fake_rates.append(fake_rate)
    lrs.append(current_lr)

    print(f"Epoch {epoch:02d}: loss={loss:.4f}  val_auc={auc:.6f}  "
          f"fake_rate@99%={fake_rate*100:.4f}%  "
          f"rolling_fr={rolling_fr*100:.4f}%  lr={current_lr:.2e}")

    # Checkpoint on val AUC
    if auc > best_val_auc + min_delta:
        best_val_auc = auc
        no_improve   = 0
        torch.save(model.state_dict(),
                   os.path.join(args.output_dir, "model_best.pt"))
        print(f"  -> new best (val_auc={auc:.6f})")
    else:
        no_improve += 1
        if no_improve >= patience:
            print(f"Early stopping at epoch {epoch} "
                  f"(no AUC improvement for {patience} epochs)")
            break

model.load_state_dict(
    torch.load(os.path.join(args.output_dir, "model_best.pt"), weights_only=True)
)
print(f"\nLoaded best model (val_auc = {best_val_auc:.6f})")

# ---------------------------
# Final eval
# ---------------------------
auc, y_true, y_pred = evaluate(model, test_loader)
print(f"Test AUC: {auc:.4f}")

real_scores_test = y_pred[y_true == 1]
fake_scores_test = y_pred[y_true == 0]

targets = [0.90, 0.95, 0.96, 0.97, 0.98, 0.99, 0.999]
wp_thresholds = {t: float(np.quantile(real_scores_test, 1 - t)) for t in targets}

thr_90  = wp_thresholds[0.90]
thr_95  = wp_thresholds[0.95]
thr_99  = wp_thresholds[0.99]
thr_999 = wp_thresholds[0.999]

WP_LINES = [
    (thr_90,  STYLE["threshold_90"],  "90%"),
    (thr_95,  STYLE["threshold_95"],  "95%"),
    (thr_99,  STYLE["threshold_99"],  "99%"),
    (thr_999, STYLE["threshold_999"], "99.9%"),
]

# Collect raw features for diagnostic plots
print("Collecting test features for diagnostics...")
test_feats_raw, test_ys = [], []
for x, y in DataLoader(Pt2Dataset(test_files, TEST), batch_size=None):
    test_feats_raw.append(x.numpy())
    test_ys.append(y.numpy())
test_feats_raw = np.concatenate(test_feats_raw)
test_ys        = np.concatenate(test_ys)

# ---------------------------
# Plots
# ---------------------------

# 1. ROC
fpr, tpr, _ = roc_curve(y_true, y_pred)
fig, ax = plt.subplots(figsize=(6, 5))
ax.plot(fpr, tpr, color=STYLE["accent"], lw=2, label=f"AUC = {auc:.4f}")
ax.plot([0,1],[0,1], "--", color=STYLE["neutral"], lw=1, alpha=0.5)
ax.fill_between(fpr, tpr, alpha=0.08, color=STYLE["accent"])
style_ax(ax, title="ROC curve", xlabel="False positive rate", ylabel="True positive rate")
ax.legend(fontsize=10)
savefig("roc.png")

# 2. Efficiency vs fake rate (log x)
fig, ax = plt.subplots(figsize=(7, 4))
ax.plot(fpr, tpr, color=STYLE["accent"], lw=2)
ax.set_xscale("log")
for thr, col, lbl in WP_LINES:
    ax.axhline(float(lbl.rstrip("%"))/100, ls="--", lw=1, color=col, label=lbl)
style_ax(ax, title="Efficiency vs fake rate", xlabel="Fake rate (log)", ylabel="Real efficiency")
ax.legend(fontsize=9)
savefig("eff_vs_fake.png")

# 3. Score distributions
for log_y, suffix in [(False, ""), (True, "_log")]:
    fig, ax = plt.subplots(figsize=(7, 4))
    bins = np.linspace(0, 1, 100)
    ax.hist(y_pred[y_true==1], bins=bins, alpha=0.6, density=True,
            color=STYLE["real_color"], label="real", edgecolor="none")
    ax.hist(y_pred[y_true==0], bins=bins, alpha=0.6, density=True,
            color=STYLE["fake_color"], label="fake", edgecolor="none")
    for thr, col, lbl in WP_LINES:
        ax.axvline(thr, color=col, lw=1.2, ls="--", label=f"thr@{lbl}")
    if log_y: ax.set_yscale("log")
    style_ax(ax, title="Score distribution", xlabel="Model score", ylabel="Density")
    ax.legend(fontsize=8)
    savefig(f"score_dist{suffix}.png")

# 4. Score overlap zoom
fig, axes = plt.subplots(1, 2, figsize=(12, 4))
zoom_bins = np.linspace(max(0, thr_99 - 0.15), min(1, thr_99 + 0.15), 80)
for ax, yscale, title in zip(axes, ["linear","log"],
                              ["Score near 99% threshold",
                               "Score near 99% threshold (log y)"]):
    ax.hist(y_pred[y_true==1], bins=zoom_bins, alpha=0.65, density=True,
            color=STYLE["real_color"], label="real", edgecolor="none")
    ax.hist(y_pred[y_true==0], bins=zoom_bins, alpha=0.65, density=True,
            color=STYLE["fake_color"], label="fake", edgecolor="none")
    ax.axvline(thr_99, color=STYLE["threshold_99"], lw=1.5, ls="--",
               label=f"thr={thr_99:.3f}")
    ax.set_yscale(yscale)
    style_ax(ax, title=title, xlabel="Score", ylabel="Density")
    ax.legend(fontsize=9)
fig.tight_layout()
savefig("score_overlap_zoom.png")

# 5. PR curve
precision, recall, _ = precision_recall_curve(y_true, y_pred)
fig, ax = plt.subplots(figsize=(6, 4))
ax.plot(recall, precision, color=STYLE["accent"], lw=2)
ax.fill_between(recall, precision, alpha=0.08, color=STYLE["accent"])
style_ax(ax, title="Precision-Recall curve", xlabel="Recall", ylabel="Precision")
savefig("pr.png")

# 6. Calibration
prob_true, prob_pred = calibration_curve(y_true, y_pred, n_bins=20)
fig, ax = plt.subplots(figsize=(5, 5))
ax.plot(prob_pred, prob_true, "o-", color=STYLE["accent"], lw=2, ms=5, label="model")
ax.plot([0,1],[0,1],"--", color=STYLE["neutral"], lw=1, label="perfect")
style_ax(ax, title="Calibration", xlabel="Mean predicted score",
         ylabel="Fraction of positives")
ax.legend(fontsize=9)
savefig("calibration.png")

# 7. Training curves
fig = plt.figure(figsize=(14, 4))
gs  = gridspec.GridSpec(1, 3, figure=fig, wspace=0.35)

ax0 = fig.add_subplot(gs[0])
ax0.plot(train_losses, color=STYLE["accent"], lw=2)
style_ax(ax0, title="Training loss", xlabel="Epoch", ylabel="Loss")

ax1 = fig.add_subplot(gs[1])
ax1.plot(val_aucs, color=STYLE["real_color"], lw=2)
ax1.set_ylim(max(0.99, min(val_aucs) - 0.0005), 1.0)
style_ax(ax1, title="Val AUC", xlabel="Epoch", ylabel="AUC")

ax2 = fig.add_subplot(gs[2])
fr_pct  = [r * 100 for r in val_fake_rates]
roll_fr = []
buf = deque(maxlen=3)
for v in fr_pct:
    buf.append(v)
    roll_fr.append(float(np.mean(buf)))
best_ep = int(np.argmin(fr_pct))
ax2.plot(fr_pct,  color=STYLE["fake_color"], lw=1.5, alpha=0.5, label="per-epoch")
ax2.plot(roll_fr, color=STYLE["fake_color"], lw=2,   label="3-ep rolling avg")
ax2.scatter([best_ep], [fr_pct[best_ep]], color=STYLE["fake_color"], s=60, zorder=5,
            label=f"best={fr_pct[best_ep]:.3f}% @ ep{best_ep}")
ax2.set_yscale("log")
style_ax(ax2, title="Val fake rate @ 99% eff", xlabel="Epoch", ylabel="Fake rate (%)")
ax2.legend(fontsize=9)

ax0r = ax0.twinx()
ax0r.plot(lrs, color=STYLE["neutral"], lw=1, ls="--", alpha=0.6)
ax0r.set_ylabel("LR", fontsize=9, color=STYLE["neutral"])
ax0r.tick_params(labelsize=8, colors=STYLE["neutral"])
ax0r.spines[["top","right"]].set_alpha(0.3)
savefig("training_curves.png")

# 8. Feature distributions
print("Plotting feature distributions...")
ncols = 4
nrows = (len(FEATURE_NAMES) + ncols - 1) // ncols
fig, axes = plt.subplots(nrows, ncols, figsize=(16, nrows * 3))
axes = axes.flatten()

for i, fname in enumerate(FEATURE_NAMES):
    ax   = axes[i]
    vals = test_feats_raw[:, i]
    lo, hi = np.percentile(vals, 1), np.percentile(vals, 99)
    bins = np.linspace(lo, hi, 60)
    ax.hist(vals[test_ys==1], bins=bins, alpha=0.6, density=True,
            color=STYLE["real_color"], label="real", edgecolor="none")
    ax.hist(vals[test_ys==0], bins=bins, alpha=0.6, density=True,
            color=STYLE["fake_color"], label="fake", edgecolor="none")
    ax.set_title(fname, fontsize=9)
    ax.set_yscale("log")
    ax.tick_params(labelsize=7)
    ax.set_facecolor(STYLE["bg"])
    ax.spines[["top","right"]].set_visible(False)
    ax.grid(True, alpha=0.2, linewidth=0.5)

for j in range(len(FEATURE_NAMES), len(axes)):
    axes[j].set_visible(False)

handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc="lower right", fontsize=10)
fig.tight_layout()
savefig("feature_distributions.png")

# 9. Score vs track quality
print("Plotting score vs track quality...")
real_mask = test_ys == 1
real_pred = y_pred[real_mask]
real_feat = test_feats_raw[real_mask]

quality_vars = {
    "pls_nhit":    FEATURE_NAMES.index("pls_nhit"),
    "ls_pt":       FEATURE_NAMES.index("ls_pt"),
    "pt2_delta_R": FEATURE_NAMES.index("pt2_delta_R"),
}

for vname, vidx in quality_vars.items():
    fig, ax = plt.subplots(figsize=(7, 4))
    sc = ax.scatter(real_feat[:, vidx], real_pred,
                    c=real_pred, cmap="RdYlGn",
                    alpha=0.25, s=2, vmin=0, vmax=1, rasterized=True)
    for thr, col, lbl in WP_LINES:
        ax.axhline(thr, color=col, lw=1.2, ls="--", label=lbl)
    style_ax(ax, title=f"Real track score vs {vname}",
             xlabel=vname, ylabel="Model score")
    ax.legend(fontsize=8)
    plt.colorbar(sc, ax=ax, label="Score")
    savefig(f"score_vs_{vname}.png")

# 10. Working point scan
wp_effs            = np.linspace(0.80, 0.9999, 200)
wp_fake_rates_scan = []
for eff in wp_effs:
    thr = np.quantile(real_scores_test, 1 - eff)
    wp_fake_rates_scan.append((fake_scores_test > thr).mean())

fig, ax = plt.subplots(figsize=(7, 4))
ax.semilogy(wp_effs * 100, wp_fake_rates_scan, color=STYLE["accent"], lw=2)
for thr, col, lbl in WP_LINES:
    ax.axvline(float(lbl.rstrip("%")), color=col, lw=1.2, ls="--", label=lbl)
style_ax(ax, title="Fake rate vs real efficiency (working point scan)",
         xlabel="Real efficiency (%)", ylabel="Fake rate")
ax.legend(fontsize=9)
ax.grid(True, which="both", alpha=0.3, linewidth=0.5)
savefig("working_point_scan.png")

# 11. Score class comparison
fig, axes = plt.subplots(1, 2, figsize=(12, 4))
bins = np.linspace(0, 1, 120)
for ax, yscale, suffix in zip(axes, ["linear","log"], ["", " (log y)"]):
    ax.hist(y_pred[y_true==1], bins=bins, alpha=0.65, density=True,
            color=STYLE["real_color"],
            label=f"real  n={int((y_true==1).sum()):,}", edgecolor="none")
    ax.hist(y_pred[y_true==0], bins=bins, alpha=0.65, density=True,
            color=STYLE["fake_color"],
            label=f"fake  n={int((y_true==0).sum()):,}", edgecolor="none")
    for thr, col, lbl in WP_LINES:
        ax.axvline(thr, color=col, lw=1.0, ls="--")
    ax.set_yscale(yscale)
    style_ax(ax, title=f"Score distribution — real vs fake{suffix}",
             xlabel="Model score", ylabel="Density")
    ax.legend(fontsize=9)
fig.tight_layout()
savefig("score_class_comparison.png")

# 12. Fake score CDF
fig, ax = plt.subplots(figsize=(7, 4))
xs = np.sort(fake_scores_test)
ys_cdf = np.arange(1, len(xs)+1) / len(xs)
ax.plot(xs, ys_cdf, color=STYLE["fake_color"], lw=2, label="Fake CDF")
ax.set_xscale("log")
for thr, col, lbl in WP_LINES:
    fr_val = (fake_scores_test > thr).mean()
    ax.axvline(thr, color=col, lw=1.2, ls="--",
               label=f"{lbl}: FR={fr_val*100:.3f}%")
    ax.axhline(1 - fr_val, color=col, lw=0.8, ls=":", alpha=0.5)
style_ax(ax, title="Fake score CDF (log x)",
         xlabel="Score threshold (log)", ylabel="Fraction of fakes passing")
ax.legend(fontsize=9)
savefig("fake_cdf.png")

# ---------------------------
# SHAP
# ---------------------------
if not args.skip_shap:
    import shap
    print("Computing SHAP values (use --skip_shap to bypass)...")

    bg_idx  = np.random.choice(len(test_feats_raw),
                               size=min(1000, len(test_feats_raw)), replace=False)
    bg_data = torch.tensor(normalize(test_feats_raw[bg_idx]),
                           dtype=torch.float32).to(DEVICE)

    r_idx = np.where(test_ys == 1)[0]
    f_idx = np.where(test_ys == 0)[0]
    explain_idx = np.concatenate([
        np.random.choice(r_idx, min(2000, len(r_idx)), replace=False),
        np.random.choice(f_idx, min(2000, len(f_idx)), replace=False),
    ])
    explain_data = torch.tensor(
        normalize(test_feats_raw[explain_idx]), dtype=torch.float32
    ).to(DEVICE)
    explain_labels = test_ys[explain_idx]
    explain_np     = explain_data.cpu().numpy()

    model.eval()
    shim = ShapShim(model).to(DEVICE)
    shim.eval()

    explainer   = shap.DeepExplainer(shim, bg_data)
    shap_values = explainer.shap_values(explain_data)

    if isinstance(shap_values, list):
        shap_arr = np.array(shap_values[0])
    else:
        shap_arr = np.array(shap_values)

    while shap_arr.ndim > 2 and shap_arr.shape[-1] == 1:
        shap_arr = shap_arr.squeeze(-1)

    assert shap_arr.ndim == 2 and shap_arr.shape == (len(explain_idx), len(FEATURE_NAMES)), \
        f"Unexpected shap_arr shape: {shap_arr.shape}"

    # 13. SHAP global importance bar chart
    mean_abs = np.abs(shap_arr).mean(axis=0)
    order    = np.argsort(mean_abs)

    fig, ax = plt.subplots(figsize=(8, 6))
    colors = [STYLE["real_color"] if v > np.median(mean_abs) else STYLE["neutral"]
              for v in mean_abs[order]]
    ax.barh([FEATURE_NAMES[i] for i in order], mean_abs[order],
            color=colors, edgecolor="none")
    style_ax(ax, title="Feature importance (mean |SHAP|)",
             xlabel="Mean |SHAP value|")
    fig.tight_layout()
    savefig("shap_importance.png")

    # 14. SHAP summary (all samples)
    plt.figure(figsize=(8, 6))
    shap.summary_plot(shap_arr, explain_np,
                      feature_names=FEATURE_NAMES, show=False, plot_size=None)
    plt.tight_layout()
    savefig("shap_summary.png")

    # 15. SHAP beeswarm — real tracks only
    real_mask_explain = explain_labels == 1
    if real_mask_explain.sum() > 10:
        plt.figure(figsize=(8, 6))
        shap.summary_plot(shap_arr[real_mask_explain],
                          explain_np[real_mask_explain],
                          feature_names=FEATURE_NAMES, show=False,
                          plot_size=None, plot_type="dot")
        plt.tight_layout()
        savefig("shap_summary_real.png")

    # 16. SHAP beeswarm — fake tracks only
    fake_mask_explain = explain_labels == 0
    if fake_mask_explain.sum() > 10:
        plt.figure(figsize=(8, 6))
        shap.summary_plot(shap_arr[fake_mask_explain],
                          explain_np[fake_mask_explain],
                          feature_names=FEATURE_NAMES, show=False,
                          plot_size=None, plot_type="dot")
        plt.tight_layout()
        savefig("shap_summary_fake.png")

    np.save(os.path.join(args.output_dir, "shap_values.npy"), shap_arr)
    print("SHAP done.")

# ---------------------------
# Working points + metrics.json
# ---------------------------
results = {}
print("\n=== Working Points ===")
for t in targets:
    threshold      = wp_thresholds[t]
    fake_eff       = float((fake_scores_test > threshold).mean())
    fake_count_est = fake_eff * 1.2e9
    print(f"{t*100:.1f}%  threshold={threshold:.6f}  "
          f"fake={fake_eff*100:.6f}%  (~{fake_count_est:,.0f} abs fakes)")
    results[f"{t*100:.1f}"] = {
        "real_efficiency":       t,
        "threshold":             float(threshold),
        "fake_rate":             fake_eff,
        "fake_rate_pct":         round(fake_eff * 100, 6),
        "est_abs_fakes_1p2B":    round(fake_count_est, 0),
    }

with open(os.path.join(args.output_dir, "metrics.json"), "w") as f:
    json.dump(results, f, indent=2)
print(f"Saved metrics.json with {len(results)} working points.")

# ---------------------------
# Save models
# ---------------------------
torch.save(model.state_dict(), os.path.join(args.output_dir, "model_final.pt"))
print("\nDone.  model_best.pt = best val AUC  |  model_final.pt = last epoch")

model.eval()
dummy = torch.randn(1, input_dim).to(DEVICE)
torch.onnx.export(
    model, dummy,
    os.path.join(args.output_dir, "model.onnx"),
    input_names=["features"],
    output_names=["score"],
    # Variable batch size, so src/pt2_scorer.cc can score many pT2s per call
    dynamic_axes={"features": {0: "batch"}, "score": {0: "batch"}},
    opset_version=13,
    do_constant_folding=True
)
print("ONNX model exported.")