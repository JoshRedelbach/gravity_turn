""" ===============================================
                GRAVITY TURN SOLVER
=============================================== """

import numpy as np
from scipy.optimize import bisect
from single_run import run

def kick_angle_objective(kick_angle, throttle):
    
    time, data = run(throttle, kick_angle)
    last_gamma = data[3,-1]
    print("Initial conditions: ", data[:,0])
    print("Final conditions", data[:,-1])
    print("Kick angle", np.rad2deg(kick_angle))
    print("Last gamma: ", np.rad2deg(last_gamma))
    print("Throttle: ", throttle)
    return last_gamma

#def second_throttle_objective():

def find_initial_kick_angle(throttle):
    return bisect(kick_angle_objective, -np.deg2rad(40), -np.deg2rad(5), xtol=1e-1, args=(throttle,))

alpha = find_initial_kick_angle(1.1)
print(alpha)