# script to run a StellaLuna ZAMS interior structure model

# this script is prepared to reproduce the interior structure
# of a 1.33 M_sun star, but can be altered to produce convective-interior,
# radiative envelope stars (e.g. M>1.1ish M_sun)

# imports
import numpy as np
from scipy.optimize import least_squares

# constants
import constants as c
from scipy.constants import N_A

# StellaLuna imports
import interpolate
import energy
import density
from guesses import *
from ode import ode
from shootf import shooter, solver

# guesses
# surface
M_star = 1.33*c.Ms
L_star_starting = ((M_star/c.Ms)**(3.5))*c.Ls # eq. 1.88 HKT
R_star_starting = ((M_star/c.Ms)**(0.75))*c.Rs # eq. 1.87 HKT
# core
Pc_starting = (3/(8*np.pi))*(c.G*(M_star)**2)/(R_star_starting)**4 # constant density sphere, lower limit!
Tc_starting = (1/2)*((4/(3+5*0.7))/(N_A*c.k))*(c.G*M_star)/(R_star_starting) # constant density sphere

print('Starting Guess:')
print('L/Lsun', L_star_starting/c.Ls)
print('logPc', np.log10(Pc_starting))
print('R/Rsun', R_star_starting/c.Rs)
print('logTc', np.log10(Tc_starting))

# initial guess vector
vec = np.array([L_star_starting, Pc_starting, R_star_starting, Tc_starting])
# shootf args
# Star mass, Shooting point in fraction of Mass, number of saved points, interior starting point, exterior starting point, do multiprocessing?
args = (M_star, 0.25, int(1e5), 1e-12, 0.9999, True)
# set limits for the minimizer
bounds = (np.array([1e-1*c.Ls, Pc_starting, 1e-1*c.Rs, Tc_starting]),
          np.array([1e6*c.Ls, Pc_starting*1e3, 1e3*c.Rs, Tc_starting*1e2]))

# run least_squares minimizer to converge model
# this is effectively the "newton" solver here
# this minimizes the difference between two shootf runs
final = least_squares(shooter, vec, args=args, bounds=bounds,
                      method='dogbox', loss='arctan',
                      gtol=None,
                      xtol=None,
                      ftol=1e-6,
                      x_scale='jac',
                     )

print(final)

if np.sum(final.active_mask**2) != 0:
    print('something ran up against a bound')

# run solution and create densely sampled results table
solution = solver(final.x, M_star=args[0], M_fit=args[1], n=1e6, in_factor=args[3], out_factor=args[4], multithread=args[5])

L_star, Pc, R_star, Tc = final.x
print('Converged starting values:')
print('L/Lsun', L_star/c.Ls)
print('logPc', np.log10(Pc))
print('R/Rsun', R_star/c.Rs)
print('logTc', np.log10(Tc))

np.savetxt('converged_start_{}.txt'.format(M_star/c.Ms), final.x)

# what is an appropriate P_factor to speed up convergence?
print('ratio between constant density Pc and converged solution',Pc/Pc_starting)

print('ratio between constant density Tc and converged solution',Tc/Tc_starting)

print('Radius is ',R_star/R_star_starting, 'of homology guess')
print('Luminosity is ',L_star/L_star_starting, 'of homology guess')

# central density vs avg density
converged_concentration = solution[5].max()/(4*np.pi*M_star/(3*R_star**3))

print("rho_c/rho_avg is", str(round(converged_concentration,2)))

# save dense results table to disk
with open('converged_interior_{}.npy'.format(M_star/c.Ms), 'wb') as f:
    np.save(f, solution)
