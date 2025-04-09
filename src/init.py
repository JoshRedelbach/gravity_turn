""" ===============================================
                  INITIALIZATION
=============================================== """

import numpy as np

# -------------- Select Rocket ----------------
LV = 'falcon9'                  # Falcon 9

# -------------- Select Simulation Type ----------------
SYM_TYPE = 4

"""
    (1 -> Single Run)
    (2 -> Single Run Full)
    (3 -> Optimization Direct Orbit Injection No Coast)
    4 -> Optimization Orbit Injection With Coasting Single Burn
    5 -> Optimization Orbit Injection With Coasting Double Burn
"""

# ===================================================
# General parameters
# ===================================================

# -------------- Gravity Turn --------------
TIME_TO_START_KICK = 7.5                        # time to start gravity turn; [s]
DURATION_INITIAL_KICK = 45.                     # duration of gravity turn; [s]

# -------------- Desired Orbit --------------
ALT_DESIRED = 500e3                             # altitude of desired orbit; [m]

# -------------- Optimization --------------
ALPHA_LOWEST = -np.deg2rad(4.)                  # lowest possible kick angle to be tested; [rad]
ALPHA_HIGHEST = -np.deg2rad(3.)                 # highest possible kick angle to be tested; [rad]
ALPHA_INITIAL_GUESS = - np.deg2rad(10.5)        # initial guess for kick angle; [rad]
ALT_NO_ATMOSPHERE = 65e3                        # altitude where atmosphere can be neglected; [m]
MAX_ACCEPTED_BURN_TIME = 15.                    # maximum accepted burn time of delta-v; [s]

OPTIMIZATION_METHOD = 2                         # optimization method
"""
    1 -> Differential Evolution
    2 -> Brute Force
"""

# ===================================================
# Single Run specific parameters
# ===================================================
SS_THROTTLE = 1.0                               # Second Stage throttle 
INITIAL_KICK_ANGLE = - np.deg2rad(1.5)          # Initial kick angle [rad]


# ===================================================
# FOR SIMULATION
# ===================================================
TIME_STEP = 0.01                              # step size for integration; [s]
DURATION_AFTER_SIMULATION = 1000.               # duration of simulation after reaching desired orbit; [s]


# ===================================================
# FOR DEBUGGING
# ===================================================
INTERRUPTS_PRINT = False
EVENTS_PRINT = False