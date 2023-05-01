# imports
import numpy as np
import pandas as pd

# plotting
import matplotlib.pyplot as plt
import seaborn as sb
sb.set_context("talk")
plt.style.use('dark_background')
# plt.style.use('default')
plt.rcParams['font.family'] = 'monospace'   # Fonts
plt.rcParams['font.monospace'] = 'DejaVu Sans Mono'

# choose X, Y, Z
# begin with X = 0.7, Z=0.02

opacity_folder = './opacities/'
# from https://www.wichita.edu/academics/fairmount_college_of_liberal_arts_and_sciences/physics/Research/opacity.php
low_T_opacity_path = opacity_folder+'A09photo.7.02.tron'
# from https://opalopacity.llnl.gov/
high_T_opacity_path = opacity_folder+'OPAL_GN93.7.02.txt'

# these opacity tables aren't that large, so we don't lose time loading them using pandas
# we get to maintain a nicer format and can pass them to quicker arrays when we need to compute with them
colw = list(np.ones(19, dtype=int)*7)
colw.insert(0, 6)
low_T_rosseland = pd.read_fwf(low_T_opacity_path, header=2, widths=colw, index_col=0)

high_T_rosseland = pd.read_csv(high_T_opacity_path, index_col=0).iloc[::-1] # reverse order to highest T at top

# set all on same log_R grid
log_R = high_T_rosseland.columns
low_T_rosseland.columns = log_R

# concatenate low and high T
rosseland = pd.concat([high_T_rosseland, low_T_rosseland], axis=0)
# drop duplicate rows where the two grids overlap, keeping the lower T rows
rosseland = rosseland[~rosseland.index.duplicated(keep='last')]

print('opacities stitched!')

rosseland.index.name = 'logT'

# save the result to the opacities folder
rosseland.to_csv(opacity_folder+'combined_OPAL_F05_X0.7_Z0.02_logT8.7-2.7.txt')
print('saved opacities to  '+opacity_folder+'combined_OPAL_F05_X0.7_Z0.02_logT8.7-2.7.txt')

# recreate figure
log_T = rosseland.index
log_R = rosseland.columns

cmap = plt.cm.magma

fig = plt.subplots(1,1,figsize=(9, 6))
for i,ind in enumerate(log_R):
    plt.plot(log_T, rosseland[ind], color=cmap(i/len(log_R)))
sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=-8, vmax=1))
plt.colorbar(sm, ax=plt.gca(), label=r'log(R)')
plt.suptitle('Roseland Mean Opacities for GN93hz mixture + F05 mixture \n(X=0.7, Y=0.28, Z=0.02)', fontsize=16)
plt.xlabel('Temperature (log K)')
plt.ylabel(r'log($\kappa_{rad}$) [cm$^{2}$/g$^{-1}$]')
plt.savefig('./figures/extended_opacity_alt.png', dpi=300, bbox_inches='tight', transparent=True)

print('saved figure to extended_opacity.png')
