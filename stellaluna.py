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
L_star = (M_star/c.Ms)**(3.5)*c.Ls # eq. 1.88 HKT
R_star = (M_star/c.Ms)**(0.75)*c.Rs # eq. 1.87 HKT
# core
Pc = (3/(8*np.pi))*(c.G*(M_star)**2)/(R_star)**4 # constant density sphere, lower limit!
Tc = (1/2)*((4/(3+5*0.7))/(N_A*c.k))*(c.G*M_star)/(R_star) # constant density sphere
P_factor = 1e1 # guess to inflate constant density sphere pressure
Pc *= P_factor

# initial guess vector
vec = np.array([L_star, Pc, R_star, Tc])
# shootf args
# Star mass, Shooting point in fraction of Mass, number of saved points, interior starting point, exterior starting point, do multiprocessing?
args = (M_star, 0.33, int(1e5), 1e-12, 0.9999, True)
# set limits for the minimizer
bounds = ([L_star*1e-1,Pc/P_factor,R_star*1e-1,Tc],
          [L_star*1e1, Pc/P_factor*1e4, R_star*1e1, Tc*1e3])

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

# what is an appropriate P_factor to speed up convergence?
print('ratio between constant density Pc and converged solution',solution[2].max()/(Pc/P_factor))

print('ratio between constant density Tc and converged solution',solution[4].max()/Tc)

L_star, Pc, R_star, Tc = final.x

print('Radius is ',solution[3].max()/((M_star/c.Ms)**(0.75)*c.Rs), 'of homology guess')
print('Luminosity is ',solution[1].max()/((M_star/c.Ms)**(3.5)*c.Ls), 'of homology guess')

# central density vs avg density
converged_concentration = solution[5].max()/(4*np.pi*M_star/(3*R_star**3))

print("rho_c/rho_avg is", str(round(converged_concentration,2)))

# save dense results table to disk
with open('converged_interior_{}.npy'.format(M_star/c.Ms), 'wb') as f:
    np.save(f, solution)
