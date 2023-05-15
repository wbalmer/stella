# imports
import numpy as np
from scipy import optimize
# stella imports
import interpolate
import energy
import density
# constants
import constants as c
from scipy.constants import N_A


def load_inner(Tc, Pc, factor=1e-12):
    """
    Returns reasonable guess for an integration
    starting point exterior to the absolute center
    """

    rho_c = density.density(Pc,Tc, X=0.7) # calculate the density of the core

    m = factor*c.Ms # set start just outside center in M_r
    epsilon = energy.energy_gen(rho_c,Tc) # determine energy generation from pp and CNO
    l = epsilon*m # calculate luminosity at core

    # KWW eq. 11.3
    r = (3/(4*np.pi*rho_c))**(1/3) * m**(1/3) # radius starting point

    # KWW eq. 11.6
    P = Pc - (3*c.G/(8*np.pi))*((4*np.pi*rho_c/3)**(4/3))*(m**(2/3)) # pressure just outside core

    # calculate temperature gradient
    del_rad = energy.del_rad(m, l, P, rho_c, Tc)

    # determine whether the core is convective or radiative
    if del_rad > energy.del_ad:
        # if convective, use KWW eq. 11.9 b
        # see also HKT eq. 7.110
        lnT = np.log(Tc) - (np.pi/6)**(1/3)*c.G*(energy.del_ad*rho_c**(4/3))/Pc*m**(2/3)
        T = np.exp(lnT)
    else:
        # if radiative, use KWW eq. 11.9 a
        # see also HKT eq. 7.109
        kappa_c = interpolate.interp_k(rho_c,Tc)
        T = (Tc**4 - (1/(2*c.a*c.c)*(3/(4*np.pi))**(2/3)*kappa_c*epsilon*rho_c**(4/3)*m**(2/3)))**(1/4)

    # return guess array
    return np.array([l, P, r, T])


def load_outer(M_star, L_star, R_star, factor=0.9999, X=0.7):
    """
    Returns reasonable guess for an integration
    starting point interior to the photosphere.
    """

    # set mu based on hydrogen fraction
    mu = 4/(3+5*X)
    # calculate surface gravity from mass, radius
    surface_g = c.G*M_star/(R_star**2)
    # calculate Teff from luminosity, radius
    Teff = (L_star/(4*np.pi*c.sb*R_star**2))**(1/4)

    # minimize the difference between rho based on opacity and rho based on equation of state
    # we need to find this rho so we can determine the photospheric opacity
    # and therefore the surface pressure, which is our final guess component
    def min_rho(rho):
        kappa = interpolate.interp_k(rho,Teff)
        # eq. 4.48 in HKT
        opacity_pressure = (2/3) * (surface_g / kappa) #* (1 + (kappa*L_star/(4*np.pi*c.c*c.G*M_star)))
        gas_rad_pressure = (1/3)*c.a*Teff**4 + (rho * N_A*c.k*Teff/mu)
        diff = 1 - opacity_pressure/gas_rad_pressure
        return np.abs(diff**2)
    # minimize this difference
    rho_sol = optimize.minimize(min_rho, 1e-8, args=(), method='Nelder-Mead', bounds=[(1e-13,1e-5)])
    # determine whether the minimizer was successful
    if rho_sol.success:
        rho = rho_sol.x[0]
    else:
        # return a nan if we interpolate off the grid or have some weird negative value
        print('there\'s no rho for this Teff, log(g)')
        rho = np.nan
    # calculate surface opacity and pressure
    kappa = interpolate.interp_k(rho,Teff)
    P = 2*surface_g/(3*kappa) #* (1 + (kappa*L_star/(4*np.pi*c.c*c.G*M_star)))
    # return guess array
    return np.array([L_star, P, R_star, Teff])
