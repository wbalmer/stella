import numpy as np
import pandas as pd
from scipy.interpolate import interpn

rosseland = pd.read_csv('./opacities/combined_OPAL_F05_X0.7_Z0.02_logT8.7-2.7.txt', index_col=0)
# make the columns floats not strings (useful later)
rosseland.columns = rosseland.columns.astype(float)

def interp_k(rho_i, T_i, method='linear'):
    '''
    A function to interpolate across the Rosseland 
    opacity grid using scipy.interpolate.interpn
    '''
    # calc logR (R=rho/T_6^3)
    logR_i = np.log10(rho_i/(T_i/1e6)**3)
    # calc logT
    logT_i = np.log10(T_i)
    # setup the grid
    points = (rosseland.index.values,rosseland.columns.values)
    values = rosseland.values
    # arrange interpolation point in array
    new_point = np.array([logT_i,logR_i])
    # compute interpolation
    return interpn(points, values, new_point, method=method)