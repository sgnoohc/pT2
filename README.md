# pT2

This repository is used to study the potential for creating and adding pT2 objects into the CMS LST algorithm.

## pT2 Objects

pT2 objects are constructed of the following:

- **pLS**: pixel line segment from the inner tracker
- **LS**: line segment from the outer tracker

## Input

As of February 06, 2026, this repository uses the output from the LST OD-ntuple.

## Building

Clone this repo and run the following in any area on HiPerGator:
(Instructions for running on uaf will come soon)

```
source setup.sh
make clean
make -j8
```

## Running

The code runs in three stages that share a histogram file (`<output dir>/pt2_hists.root`):

```
./bin/pt2 process -k -r      # event loop: build pT2s, fill histograms (-r also writes LSTNtuple_with_pT2.root)
./bin/pt2 scan -e 90         # compute cut values at the target efficiency from the histograms
./bin/pt2 plot               # draw plots + index.html from the histograms
```

`scan` and `plot` only read the histogram file, so they can be rerun without redoing the event loop.

Options:

```
process
  -i Input File Path
  -k low pt (we usually use this)
  -n Number of events
  -r make root file adding pt2s
scan
  -e Target efficiency percent (default 90)
common
  -o Output Directory (default: output)
  -H Histogram file (default: <output dir>/pt2_hists.root)
```

## Training the NN

`process -r` also writes `<output dir>/pt2_training_data.root`: a flat tree with one row per pT2 and an `event_idx` branch.
Train on one or more of these files with the PyTorch module (in a fresh shell, not one where `setup.sh` was sourced):

```
module load pytorch/2.8.0
python pt2_ml/train_v7.py --data output/pt2_training_data.root --output_dir nn_out --skip_shap
```

- `--data` takes files, directories, or glob patterns. Rows are split into train/val/test by event (`--val_frac`, `--test_frac`, default 0.1 each).
- `--fast_dev_run` does a quick 2-epoch check on a few batches.
- SHAP plots need the `shap` package, which the module does not include; use `--skip_shap` otherwise.
- The output directory gets `model.onnx`, `mean.npy`, `std.npy` (exported with a dynamic batch size), plus plots and `metrics.json`. Use it with `./bin/pt2 process -r -N nn_out`.
