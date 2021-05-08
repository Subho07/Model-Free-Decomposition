This package performs model free 3-component scattering power decomposition for full polarimetric Synthetic Aperture Radar data.

Dey, S., Bhattacharya, A., Ratha, D., Mandal, D. and Frery, A.C., 2020. Target Characterization and Scattering Power Decomposition for Full and Compact Polarimetric SAR Data. IEEE Transactions on Geoscience and Remote Sensing.

## Python usage

from mf3cf import mf3cf_powers

ps, pd, pv = mf3cf_powers('D:/T3', 7)

path = 'D:/T3'

window_size = 7

## C usage
To compile and run with T3 matrix data in T3 folder and analysis window size 3:
```
   gcc mf3cc.c -o mf3cc.exe -lm
   ./mf3cc.exe T3 3
```
