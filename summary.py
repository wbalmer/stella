# create figures visualizing a converged StellaLuna model

# imports
import os
import numpy as np
import pandas as pd
from scipy.optimize import least_squares

# constants
import constants as c
from scipy.constants import N_A

# plotting
import matplotlib.pyplot as plt
import seaborn as sb
plt.style.use('dark_background')
# plt.style.use('default')
sb.set_context("talk")
plt.rcParams['font.family'] = 'monospace'   # Fonts
plt.rcParams['font.monospace'] = 'DejaVu Sans Mono'

plot_dir = './figures/'
if not os.path.isdir(plot_dir):
    os.makedirs(plot_dir)
    print("created folder : ", plot_dir)

# guesses
# surface
M_star = 1.33*c.Ms
lpower = 3.5
L_star_starting = (M_star/c.Ms)**(lpower)*c.Ls # eq. 1.88 HKT
rpower = 0.75
R_star_starting = (M_star/c.Ms)**(rpower)*c.Rs # eq. 1.87 HKT
# core
Pc_starting = (3/(8*np.pi))*(c.G*(M_star)**2)/(R_star_starting)**4 # constant density sphere, lower limit!
Tc_starting = (1/2)*((4/(3+5*0.7))/(N_A*c.k))*(c.G*M_star)/(R_star_starting) # constant density sphere

# load solution
solution = np.load('converged_interior_{}.npy'.format(M_star/c.Ms))
L_star, Pc, R_star, Tc = np.loadtxt('converged_start_{}.txt'.format(M_star/c.Ms))

cols = ['m', 'l', 'P', 'r', 'T', 'd', 'e', 'k', 'ad', 'rad']
data = pd.DataFrame(solution.T, columns=cols)

data['m'] /= c.Ms
data['l'] /= c.Ls
data['P'] = np.log10(data['P'])
data['r'] /= c.Rs
data['T'] = np.log10(data['T'])
data['d'] = np.log10(data['d'])
# data['e'] =
# data['k'] =
# data['ad'] =
# data['rad'] =
data['nabla'] = np.minimum(data['rad'], data['ad'])
data['trans'] = 'Radiative'
data['trans'][data['nabla']==0.4] = 'Convective'

summary = pd.concat([pd.concat([pd.concat([data.iloc[::int(1e5), :],
                                           data.iloc[1990000::int(1e4), :]]),
                                           data.iloc[1999000::int(1e3), :]]),
                                           data.iloc[1999990:1999998:3, :]]
                   ).drop_duplicates().reset_index(drop=True)

summary_text = open("latex_summary_{}.txt".format(M_star/c.Ms), "w")

summary_text.write(summary.to_latex(index=False,
                   formatters={"name": str.upper},
                   float_format="{:.3e}".format))

summary_text.close()
print('wrote summary latex table!')

summary.to_csv('summary_interior_{}.csv'.format(M_star/c.Ms), index=False)
print('wrote summary csv table!')

# compare to MESA

# R_star
stella_rstar = R_star/c.Rs
# Pc
stella_pc = np.log10(Pc)
# Tc
stella_tc = np.log10(Tc)
# L_star
stella_lstar = L_star/c.Ls
print('StellaLuna gives:')
print(stella_lstar, stella_pc, stella_rstar, stella_tc)

# mesa for 1.33 Msun
# mesa log Tc
mesa_tc = 7.200495
# mesa log Dc
mesa_rho = 10**1.819171
def pressure(rho,T, X=0.7):
    '''
    Calculate pressure assuming an EOS that is a simple mixture
    of ideal gas and radiation (e.g. HKT eq 3.104)
    '''
    # we assume X=0.7, and full ionization, so that mu ~ 2/(1+X)
    mu = 4/(3+5*X)
    return (rho*c.k*N_A*T/mu) + (c.a/3*(T**4))
mesa_pc = np.log10(pressure(mesa_rho, 10**mesa_tc))
# mesa log R
mesa_rstar = (10**(0.135255))#*c.Rs
# mesa log L
mesa_lstar = (10**(0.428695))#*c.Ls

print('MESA gives:')
print(mesa_lstar, mesa_pc, mesa_rstar, mesa_tc)

percenterr = lambda x,y : (x-y)/y * 100

pl = percenterr(stella_lstar, mesa_lstar)
pp = percenterr(stella_pc, mesa_pc)
pr = percenterr(stella_rstar,mesa_rstar)
pt = percenterr(stella_tc,mesa_tc)

print('percent error:')
print(pl, pp, pr, pt)
