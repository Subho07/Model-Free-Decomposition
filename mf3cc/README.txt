This package performs model free 3-component scattering power decomposition for compact polarimetric Synthetic Aperture Radar data.

Dey, S., Bhattacharya, A., Ratha, D., Mandal, D. and Frery, A.C., 2020. Target Characterization and Scattering Power Decomposition for Full and Compact Polarimetric SAR Data. IEEE Transactions on Geoscience and Remote Sensing.

## usage

from mf3cc import mf3cc_powers

ps, pd, pv = mf3cc_powers('D:/C2', 7, -45)

path = 'D:/C2' input a C2 folder

window_size = 7

chi_value = -45 (-45 is default value. However, you may change it as per requirement)
