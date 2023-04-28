import numpy as np
import pandas as pd
from scipy.interpolate import griddata

rosseland = pd.read_csv('./opacities/combined_OPAL_F05_X0.7_Z0.02_logT8.7-2.7.txt', index_col=0)
# rosseland = pd.read_csv('./opacities/OPAL_GN93.7.02.txt', index_col=0)
# make the columns floats not strings (useful later)
rosseland.columns = rosseland.columns.astype(float)

# setup the grid
points = []
values = []
for i in range(len(rosseland.columns.values)):
    for j in range(len(rosseland.index.values)):
        points.append([rosseland.columns.values[i], rosseland.index.values[j]])
        values.append(rosseland.values[j][i])

def interp_k(rho_i, T_i, method='linear'):
    '''
    A function to interpolate across the Rosseland
    opacity grid using scipy.interpolate.interpn
    '''
    # calc logR (R=rho/T_6^3)
    logR_i = np.log10(rho_i/(T_i/1e6)**3)
    # calc logT
    logT_i = np.log10(T_i)

    # arrange interpolation point in array
    new_point = (logT_i,logR_i)
    # compute interpolation
    interpolated_points = 10**griddata(points, values, (logR_i, logT_i), method=method)
    return interpolated_points
