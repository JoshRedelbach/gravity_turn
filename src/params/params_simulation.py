""" ===============================================
            PARAMETERS FOR SIMULATION
=============================================== """

import numpy as np

# -------------- Gravitiy Turn --------------
max_angle_of_attack = - np.deg2rad(14.8)       # max. angle of attack during pitch programm; [rad]
alt_initial_kick = 1000                     # altitude to start gravity turn; [m]
duration_initial_kick = 20.                 # duration of gravity turn; [s]


# -------------- For Integration --------------
number_of_points = 1000                     # number of data points that are calculated in scipy.solve_ivp(...)
time_step = 0.001                            # step size for integration; [s]


# -------------- Final Conditions --------------
# NOTE: TO BE IMPLEMENTED
# Desired Orbit for Final Conditions Check
alt_desired = 100e3   # altitude of desired orbit; [m]
# ...

