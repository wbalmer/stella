# imports
import numpy as np
from scipy.integrate import solve_ivp
import ray
from ray.util.multiprocessing import Pool
# stella imports
import interpolate
import energy
import density
from guesses import *
from ode import ode
# constants
import constants as c
from scipy.constants import N_A


def shooter(vec, M_star=1.67*c.Ms, M_fit=0.5,
            n=int(1e5), in_factor=1e-12, out_factor=0.9999,
            multithread=False,
            ):
    """
    a version of the "shootf" function that takes a vector of initial guesses
    (luminosity, central pressure, radius, central temperature) and shoots towards
    a solution from the interior and surface.
    set mass of star with M_star
    set start points (in fraction of enclosed mass) via in_factor and out_factor
    set fitting point (in fraction of enclosed mass) via M_fit
    Multithread the ODE integration using ray/pool can be toggled
    """

    # load in vector of initial guess variables
    L_star, Pc, R_star, Tc = vec

    # load initial guess vectors based on input variables
    inn = load_inner(Tc, Pc, factor=in_factor)
    outt = load_outer(M_star, L_star, R_star, factor=out_factor)

    # protect against bad solutions which crash the minimizer
    if np.isnan(np.sum(inn)) or np.isnan(np.sum(outt)):
        print('caught a nan in the guess')
        return np.array([-np.inf, -np.inf, -np.inf, -np.inf])

    # set up array of enclosed mass to solve across
    exiting = np.logspace(np.log10(in_factor*c.Ms), np.log10(M_fit*M_star), base = 10.0, num = int(n))
    entering = np.flipud(np.logspace(np.log10(M_fit*M_star), np.log10(M_star), base = 10.0 , num = n))

    # set up multithreading
    if multithread:
        ray.init(num_cpus=4)
        pool = Pool()
        # solve heading from core to surface
        sol_i = pool.apply(solve_ivp, [ode, (exiting[0], exiting[-1]), inn, 'RK45', exiting])
    else:
        # solve heading from core to surface
        sol_i = solve_ivp(ode, (exiting[0], exiting[-1]), inn, method='RK45', t_eval=exiting)
    # determine success of core->surface integrator
    if sol_i.status == 0:
        # report success
        print('solved inner')
    else:
        # report failure, shutdown multithread if using multithreading
        print('failed to solve interior', sol_i.message)
        if multithread:
            ray.shutdown()
        # protect minimizer against failed solutions
        return np.array([-np.inf, -np.inf, -np.inf, -np.inf])

    if multithread:
        # solve heading from surface to core
        sol_s = pool.apply(solve_ivp, [ode, (entering[0], entering[-1]), outt, 'RK45', entering])
    else:
        # solve heading from surface to core
        sol_s = solve_ivp(ode, (entering[0], entering[-1]), outt, method='RK45', t_eval=entering)
    # determine success of core->surface integrator
    if sol_s.status == 0:
        # report success
        if multithread:
            ray.shutdown()
        print('solved exterior')
    else:
        # report failure, shutdown multithread if using multithreading
        print('failed to solve exterior', sol_s.message)
        if multithread:
            ray.shutdown()
        # protect minimizer against failed solutions
        return np.array([-np.inf, -np.inf, -np.inf, -np.inf])

    # assign integrated solution to variables
    exiting_sol = sol_i.y
    entering_sol = sol_s.y

    # determine the difference at the shooting point
    dL = (exiting_sol[0,-1] - entering_sol[0,-1])/L_star
    dP = (exiting_sol[1,-1] - entering_sol[1,-1])/Pc
    dR = (exiting_sol[2,-1] - entering_sol[2,-1])/R_star
    dT = (exiting_sol[3,-1] - entering_sol[3,-1])/Tc
    # return residual array
    print(np.array([dL, dP, dR, dT]))
    return np.array([dL, dP, dR, dT])


def solver(vec_final, M_star=1.67*c.Ms, M_fit=0.5,
            n=int(1e5), in_factor=1e-12, out_factor=0.9999,
            multithread=False,
            ):
    """
    a version of the "shootf" function that solves the ode given the results of a
    minimizer and returns a final solution array
    set mass of star with M_star
    set start points (in fraction of enclosed mass) via in_factor and out_factor
    set fitting point (in fraction of enclosed mass) via M_fit
    Multithread the ODE integration using ray/pool can be toggled
    """

    L_star, Pc, R_star, Tc = vec_final

    inn = load_inner(Tc, Pc, factor=in_factor)
    outt = load_outer(M_star, L_star, R_star, factor=out_factor)

    exiting = np.logspace(np.log10(in_factor*c.Ms), np.log10(M_fit*M_star), base = 10.0, num = int(n))
    entering = np.flipud(np.linspace(M_fit*M_star, M_star, num = int(n)))

    if multithread:
        ray.init(num_cpus=4)
        pool = Pool()
        sol_i = pool.apply(solve_ivp, [ode, (exiting[0], exiting[-1]), inn, 'RK45', exiting])
        sol_s = pool.apply(solve_ivp, [ode, (entering[0], entering[-1]), outt, 'RK45', entering])
        ray.shutdown()
    else:
        sol_i = solve_ivp(ode, (exiting[0], exiting[-1]), inn, method='RK45', t_eval=exiting)
        sol_s = solve_ivp(ode, (entering[0], entering[-1]), outt, method='RK45', t_eval=entering)

    exiting_sol = sol_i.y
    entering_sol = sol_s.y

    # combine mass arrays
    mass = np.concatenate([exiting, np.flipud(entering)], axis=0)

    # add mass to final array
    solution = np.zeros((7, mass.shape[0]))
    solution[0] = mass

    # combine solution arrays
    sols = np.concatenate([exiting_sol, np.fliplr(entering_sol)], axis=1)
    solution[1:5] = sols

    # add density as 6th column
    rho = density.density(solution[2],solution[4], X=0.7)
    solution[5] = rho

    # add del_rad as 7th column
    del_rad = energy.del_rad(mass, solution[1], solution[2], rho, solution[4])
    solution[6] = del_rad

    return solution
