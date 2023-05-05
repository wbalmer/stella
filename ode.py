# imports
import numpy as np
# stella imports
import interpolate
import energy
import density
# constants
import constants as c


def ode(m,v):
    """
    The four coupled differential equations for a simple stellar model
    See KWW Ch 10 or HKT Ch 7.
    """

    # load vector of variables
    l, P, r, T = v
    # calculate density
    rho = density.density(P,T, X=0.7)
    # calculate del_rad
    del_radiative = energy.del_rad(m, l, P, rho, T)
    # determine whether the star is convective or not
    del_actual = np.minimum(del_radiative, energy.del_ad)

    # calculate the derivatives
    dldm = energy.energy_gen(rho,T) #change in luminosity with enclosed mass
    dPdm = -c.G*m/(4*np.pi*r**4) # change in pressure with enclosed mass
    drdm = 1/(4*np.pi*r**2 * rho) # mass conservation eq.
    dTdm = ((-c.G*m*T)/(4*np.pi*P*r**4))*del_actual # change in temperature with enclosed mass

    # return derivatives
    return np.array([dldm, dPdm, drdm, dTdm])
