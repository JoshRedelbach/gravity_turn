""" ===============================================
                GRAVITY TURN SOLVER
=============================================== """

import numpy as np
from scipy.optimize import bisect
from components.rocket import run

import params.constants as c
import params.params_simulation as par_sim
import components.rocket as rocket

def kick_angle_objective(kick_angle, throttle):
    time, data = run(throttle, kick_angle)
    last_gamma = data[3,-1]
    return last_gamma

def find_initial_kick_angle(throttle):
    return bisect(kick_angle_objective, -np.deg2rad(60), -np.deg2rad(0.1), xtol=1e-7, args=(throttle,), maxiter=1000)

def second_throttle_objective(throttle):
    print("----------------------------")
    print("Throttle: ", throttle)
    initial_kick_angle = find_initial_kick_angle(throttle)
    time, data = run(throttle, initial_kick_angle)
    last_radius = data[1, -1]
    
    last_gamma = data[3,-1]
    # print("Initial conditions: ", data[:,0])
    # print("Final conditions", data[:,-1])
    print("Kick angle", np.rad2deg(initial_kick_angle))
    print("Last gamma: ", np.rad2deg(last_gamma))
    print("Throttle: ", throttle)
    a, e, r_apo, r_peri = rocket.get_orbital_elements(data[1, -1], data[2, -1], data[3, -1])
    print("Perigee: ", r_peri - c.r_earth, "m")
    print("Apogee: ", r_apo - c.r_earth, "m")
    print("Eccentricity: ", e)
    print("Altitude: ", (last_radius - c.r_earth)/1000, "km")
    r_desired = c.r_earth + par_sim.alt_desired
    v_desired = np.sqrt(c.mu_earth / r_desired)
    print("delta_Velocity: ", data[2, -1] - v_desired, "m/s")
    print(last_radius - c.r_earth - par_sim.alt_desired)
    
    return last_radius - c.r_earth - par_sim.alt_desired

def find_throttle():
    return bisect(second_throttle_objective, 2, 0.5, xtol=1e-6, maxiter=500)

find_throttle()
