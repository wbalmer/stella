import numpy as np

def f11(rho,T):
    '''
    Calculate weak screening for a temperature and density,
    using the relation quoted in KWW Ch 18.4, eqs 18.56-57
    '''
    T7 = T/1e7
    term1 = 5.92*(10**(-3))
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
    return 1 + 3.82*T9 + 1.51*T9**2 + 0.144*T9**3 - 0.0114*T9**4

def g141(T):
    '''
    Calculate the gaunt factor for the CNO-cycle based on eq 18.65 in KWW CH 18.5
    '''
    T9 = T/1e9
    return 1 - 2.00*T9 + 3.41*T9**2 - 2.43*T9**3

def epsilon_pp(rho,T): 
    T9 = T/1e9
    X1 = 0.7 # assuming at the beginning based on our opacities
    psi = 1 # varies between 1 ~ 1.5 depending on T and X_He. couldn't find parameterized version so setting to 1
    f_weak = f11(rho,T)
    g_pp = g11(T)
    return 2.57*10**4 * psi * f_weak * g_pp * rho * X1**2 * T9**(-2/3) * np.exp(-3.381/(T9**(1/3)))

def epsilon_cno(rho,T):
    T9 = T/1e9
    X1 = 0.7 # assuming at the beginning based on our opacities
    X_cno = 0.01
    psi = 1 # varies between 1 ~ 1.5 depending on T and X_He. couldn't find parameterized version so setting to 1
    g_cno = g141(T)
    return 8.24*10**25*g_cno*X_cno*X1*rho*T9**(-2/3)*np.exp((-15.231*T9**(-1/3))-(T9/0.8)**2)
  

def energy_gen(rho,T):
    '''
    Calculate rate of energy generation in cgs following KWW Ch 18.5, which draws from Angulo+99
    '''

    e_pp = epsilon_pp(rho,T)
    e_cno = epsilon_cno(rho,T)
   
    return e_pp+e_cno

if __name__ == '__main__':

    x = np.linspace(6, 8, 50)
    y = np.log10(energy.energy_gen(1,10**x))
    y_pp = np.log10(energy.epsilon_pp(1,10**x))
    y_cno = np.log10(energy.epsilon_cno(1,10**x))

    plt.plot(x,y, color='k')
    plt.plot(x,y_pp, ls='--', color='k')
    plt.plot(x,y_cno, ls='--', color='k')
    plt.ylim(-10,10)

    plt.xlabel(r'$\log_{10}(T/K)$')
    plt.ylabel(r'$\log_{10}(\epsilon_H)$')

    plt.savefig('energy_generation.png', dpi=300, bbox_inches='tight')
