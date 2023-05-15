import numpy as np
import constants as c
from scipy.constants import N_A

def density(P,T, X=0.7):
    '''
    Calculate density assuming an EOS that is a simple mixture
    of ideal gas and radiation (e.g. HKT eq 3.104)
    '''
    # we assume X=0.7, and full ionization, so that mu ~ 2/(1+X)
    mu = 4/(3+5*X)
    return (P - ((c.a/3)*T**4))*mu/(N_A * c.k*T)
