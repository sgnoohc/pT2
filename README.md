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

Once correctly compiled, the code can be run as follows:

```
./bin/pt2
```

The full list of run options is as follows:

```
-i Input File Path
-o Output Directory
-k low pt (we usually use this)
-p plots
-r make root file adding pt2s
```
