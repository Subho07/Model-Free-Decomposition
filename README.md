# Model Free Decomposition techniques using SAR data

[![DOI](https://img.shields.io/badge/DOI-https%3A%2F%2Fdoi.org%2F10.1109%2FTGRS.2020.3010840-brightgreen)](https://doi.org/10.1109/TGRS.2020.3010840)

This repository contains python packages for model free decompositions using Synthetic Aperture Radar (SAR) data



## Installation

Please clone the github repository

> git clone https://github.com/Subho07/Model-Free-Decomposition.git

Go to Model-Free-Decomposition

> cd Model-Free-Decomposition

For Model free 3 component full polarimetric decomposition (MF3CF)

> cd mf3cf 

or, 

> cd Model-Free-Decomposition\mf3cf

For Model free 3 component compact polarimetric decomposition (MF3CC)

> cd mf3cc

or, 

> cd Model-Free-Decomposition\mf3cc

Install packages using

> python setup.py install

## Usage MF3CF

### MF3CF in Python
```python
>> T3_path = 'D:/T3'

>> wsi = 7 # wsi = window size

>> from mf3cf import mf3cf_powers

>> ps, pd, pv = mf3cf_powers(T3_path, wsi)
```
Please check the `T3_path` for the exported power components in `.bin` format.

### MF3CF in C
```bash
gcc mf3cf.c -o mf3cf.exe -lm
./mf3cf.exe T3
```

## Usage MF3CC

### MF3CC in Python:
```python
>> C2_path = 'D:/C2'

>> wsi = 7 # wsi = window size

>> from mf3cc import mf3cc_powers

>> ps, pd, pv = mf3cc_powers(C2_path, wsi, -45) # here chi = -45 (as default)
```
Please check the `C2_path` for the exported power components in `.bin` format.

## Reference

- S. Dey, A. Bhattacharya, D. Ratha, D. Mandal and A. C. Frery, "Target Characterization and Scattering Power Decomposition for Full and Compact Polarimetric SAR Data," in IEEE Transactions on Geoscience and Remote Sensing, doi: https://doi.org/10.1109/TGRS.2020.3010840.
