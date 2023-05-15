import numpy as np
import constants as c
import interpolate


del_ad = 0.4 # assuming ideal gas, complete ionization


def del_rad(m, l, P, rho, T):
    """
    Calculate del_rad based on 4.30 in HKT
    """
    try:
        kappa = interpolate.interp_k(rho, T)
    except ValueError:
        print(np.log10(rho/(T/1e6)**3), np.log10(T))
        raise ValueError('trying to interpolate out of bounds')
    return (3/(16*np.pi*c.a*c.c))*(P*kappa/T**4)*(l/(c.G*m))


def f11(rho,T):
    '''
    Calculate weak screening for a temperature and density,
    using the relation quoted in KWW Ch 18.4, eqs 18.56-57
    '''
    T7 = T/1e7
    term1 = 5.92e-3
    Z1Z2 = 1 # assuming pp-chain
    zeta = 1 # KWW say this is of order unity
    term2 = (zeta*rho/(T7**3))**(1/2)
    edkt = term1*Z1Z2*term2
    return np.exp(edkt)


def g11(T):
    '''
    Calculate the gaunt factor for the pp-chain based on eq 18.63 in KWW CH 18.5
    '''
    T9 = T/1e9
    return 1 + (3.82* T9) + (1.51* T9**2) + (0.144* T9**3) - (0.0114* T9**4)


def g141(T):
    '''
    Calculate the gaunt factor for the CNO-cycle based on eq 18.65 in KWW CH 18.5
    '''
    T9 = T/1e9
    return 1 - (2.00*T9) + (3.41*T9**2) - (2.43*T9**3)


def poor_mans_psi(T):
    '''
    A piecewise approximation of the pp-chain He dependency, fit by eye from
    figure 18.7 in KWW
    '''
    T7 = T / 1e7
    if T7 < 1.:
        p = 1.
    elif T7 > 3.:
        p = 1.4
    else:
        p = .2 * T7 + 0.8
    return p


def epsilon_pp(rho,T):
    '''
    energy generation for the pp-chain, from KWW Ch 18.4
    '''
    T9 = T/1e9
    X1 = 0.7 # assuming at the beginning based on our opacities
    f_weak = f11(rho,T)
    try:
        psi_T = poor_mans_psi(T)
    except ValueError:
        psi_T = np.zeros_like(T)
        for i,t in enumerate(T):
            psi_T[i] += poor_mans_psi(t)
    g_pp = g11(T)
    return 2.57e4 * psi_T * f_weak * g_pp * rho * X1**2 * T9**(-2/3) * np.exp(-3.381/(T9**(1/3)))


def epsilon_cno(rho,T,Z=0.02):
    '''
    energy generation rate for the CNO cycle, from KWW Ch 18.4
    '''
    T9 = T/1e9
    X1 = 0.7 # assuming at the beginning based on our opacities
    X_cno = (2/3)*Z
    g_cno = g141(T)
    return 8.24e25*g_cno*X_cno*X1*rho*T9**(-2/3)*np.exp((-15.231*T9**(-1/3))-(T9/0.8)**2)


def energy_gen(rho,T):
    '''
    Calculate rate of energy generation in cgs following KWW Ch 18.5, which draws from Angulo+99
    '''

    e_pp = epsilon_pp(rho,T)
    e_cno = epsilon_cno(rho,T)

    return e_pp+e_cno


if __name__ == '__main__':
    # plotting
    import matplotlib.pyplot as plt
    import seaborn as sb
    sb.set_context("talk")
    # plt.style.use('dark_background')
    plt.style.use('default')
    plt.rcParams['font.family'] = 'monospace'   # Fonts
    plt.rcParams['font.monospace'] = 'DejaVu Sans Mono'
    # this reproduces figure 18.8 in KWW, showing the contribution from each nuclear process to the total energy generation
    x = np.linspace(6, 8, 50)
    y = np.log10(energy_gen(1,10**x))
    y_pp = np.log10(epsilon_pp(1,10**x))
    y_cno = np.log10(epsilon_cno(1,10**x))

    plt.plot(x,y, color='k')
    plt.plot(x,y_pp, ls='--', color='tomato', label='pp-chain')
    plt.plot(x,y_cno, ls='--', color='cornflowerblue', label='CNO-cycle')
    plt.ylim(-10,10)

    plt.xlabel(r'$\log_{10}(T/K)$')
    plt.ylabel(r'$\log_{10}(\epsilon_H)$')
    plt.legend()

    plt.savefig('./figures/energy_generation.png', dpi=300, bbox_inches='tight', transparent=True)
