# create figures visualizing a converged StellaLuna model

# imports
import os
import numpy as np
from scipy.optimize import least_squares

# constants
import constants as c
from scipy.constants import N_A

# plotting
import matplotlib.pyplot as plt
import seaborn as sb
# plt.style.use('dark_background')
plt.style.use('default')
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
L_star_starting = (M_star/c.Ms)**(3.5)*c.Ls # eq. 1.88 HKT
R_star_starting = (M_star/c.Ms)**(0.75)*c.Rs # eq. 1.87 HKT
# core
Pc_starting = (3/(8*np.pi))*(c.G*(M_star)**2)/(R_star_starting)**4 # constant density sphere, lower limit!
Tc_starting = (1/2)*((4/(3+5*0.7))/(N_A*c.k))*(c.G*M_star)/(R_star_starting) # constant density sphere

# load solution
solution = np.load('converged_interior_{}.npy'.format(M_star/c.Ms))
L_star, Pc, R_star, Tc = np.loadtxt('converged_start_{}.txt'.format(M_star/c.Ms))

# first plot
# a figure like 9.1 in HKT but in mass space, y-values near 1
plt.figure(figsize=(9,6))

plt.plot(solution[0]/M_star, solution[5]/solution[5].max(), color='cornflowerblue', label=r'$\rho~[g/cm^3]$'.format(), ls='-.')
plt.plot(solution[0]/M_star, solution[2]/solution[2].max(), color='tomato', label=r'$P~[Ba]$', ls=':')
plt.plot(solution[0]/M_star, solution[4]/solution[4].max(), color='xkcd:forest green', label=r'$T~[K]$', ls='--')
plt.plot(solution[0]/M_star, solution[1]/solution[1].max(), color='xkcd:lavender', label=r'$L~[erg/s]$', ls='-')

plt.xlabel(r'M$_\mathrm{r}$ [M$_\star$]')
plt.ylabel(r'Normalized Quantity')
plt.legend(fontsize=13, bbox_to_anchor=(0.55,0.5))
plt.xlim(-0.1,1.1)
plt.minorticks_on()

plt.savefig(plot_dir+'run_over_mass_norm_{}.png'.format(M_star/c.Ms), dpi=300, bbox_inches='tight', transparent=True)

# second plot
# reproduce figure 9.1 in HKT but in mass space

plt.figure(figsize=(9,6))

plt.plot(solution[3]/R_star, solution[5]*1e5, color='cornflowerblue', label=r'$\log(\rho)\times10^5$ $[g/cm^3]$', ls='-.')
plt.plot(solution[3]/R_star, solution[2]*1e-8, color='tomato', label=r'$\log(P)\times10^{-8}$ $[Ba]$', ls=':')
plt.plot(solution[3]/R_star, solution[4], color='xkcd:forest green', label=r'$\log(T)$ $[K]$', ls='--')

plt.yscale('log')
plt.xlabel(r'R [R$_\star$]')
plt.ylabel(r'$\log(quantity)\times C$')
plt.legend(fontsize=10)
plt.ylim(1,1e10)
plt.xlim(0,1.1)

plt.savefig(plot_dir+'run_over_radius_HKT9-1_norm_{}.png'.format(M_star/c.Ms), dpi=300, bbox_inches='tight', transparent=True)


# third plot
# del rad vs del ad
plt.figure(figsize=(9,6))

plt.plot(solution[0]/M_star, solution[-1], color='cornflowerblue', label=r'$\nabla_{rad}$', ls=':')
plt.plot(solution[0]/M_star, np.zeros_like(solution[0])+0.4, color='tomato', label=r'$\nabla_{ad}$')

# plt.yscale('log')
plt.xlabel(r'M [M$_\star$]')
plt.ylabel(r'$\nabla$ [K/cm]')
plt.legend()
plt.ylim(0,0.6)
plt.xlim(-0.1,1.1)
plt.minorticks_on()
plt.savefig(plot_dir+'del_over_mass_norm_{}.png'.format(M_star/c.Ms), dpi=300, bbox_inches='tight', transparent=True)


# fourth plot
# central density vs avg density
converged_concentration = solution[5].max()/(4*np.pi*M_star/(3*R_star**3))

from scipy.integrate import solve_ivp

# compare to polytrope

def dtheta(xi, y, n):
    theta, psi = y
    return (psi, -np.power(theta, n) - 2 * psi / xi)

ns = np.linspace(0., 3.5, num=101) # a reasonable number of n for plotting
ts = np.logspace(np.log10(1e-15), np.log10(1000.), num=1000) # dense enough to look smooth

concentrations = np.zeros(ns.shape[0])
for i, n in enumerate(ns):
    soln = solve_ivp(dtheta, [1e-16, 1001.], [1., 0.], args=[n], t_eval=ts)
    # degree of central concentration
    concentration = 1/3 * (soln.t[np.argmin(np.abs(soln.y[0]))]/(-1*soln.y[1][np.argmin(np.abs(soln.y[0]))]))
#     print("for n={}, rho_c/rho_avg is".format(n), str(round(concentration,2)))
    concentrations[i] += concentration

x = ns
y = concentrations

points = np.array([x, y]).T.reshape(-1, 1, 2)
segments = np.concatenate([points[:-1], points[1:]], axis=1)

fig, axs = plt.subplots(1, 1, figsize=(9,5))

# Create a continuous norm to map from data points to colors
norm = plt.Normalize(x.min(), x.max())
from matplotlib.collections import LineCollection
lc = LineCollection(segments, cmap='magma', norm=norm)
# Set the values used for colormapping
lc.set_array(x)
lc.set_linewidth(2)
line = axs.add_collection(lc)
# cb = fig.colorbar(line, ax=axs)
# cb.ax.set_title('n')
plt.hlines(converged_concentration, 0, 3.5, ls='--', color='k')
plt.annotate('Model degree of \n central concentration', (2.3,6.5), fontsize=10)

plt.minorticks_on()
plt.ylim(1, 1e2)
plt.yscale('log')

plt.xlabel('Polytropic Index n')
plt.ylabel(r'$\rho_c~/~\left|\rho\right|}$')

plt.savefig(plot_dir+'central_concentration_polytropes_norm_{}.png'.format(M_star/c.Ms), dpi=300, bbox_inches='tight', transparent=True)
