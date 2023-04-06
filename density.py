import numpy as np

def density(P,T):
    '''
    Calculate density assuming an EOS that is a simple mixture 
    of ideal gas and radiation (e.g. HKT eq 3.104)
    '''
    # we assume X=0.7, and full ionization, so that mu ~ 2/(1+X)
    mu = 2/(1+0.7)
    return ((c.a/3)*T**4 - P)*mu/(N_A * c.k*T)