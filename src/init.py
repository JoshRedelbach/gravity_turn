""" ===============================================
                  INITIALIZATION
=============================================== """

import numpy as np

# -------------- Select Rocket ----------------
# LV = 'MK1'                      # Example
LV = 'falcon9'                  # Falcon 9
# LV = 'soyuz'                  # Soyuz

# -------------- Select Simulation Type ----------------

SYM_TYPE = 5

"""
    1 -> Single Run
    2 -> Single Run Full
    3 -> Optimization Direct Orbit Injection No Coast
    4 -> Optimization Orbit Injection With Coasting Single Burn
    5 -> Optimization Orbit Injection With Coasting Double Burn
"""

# ===================================================
# General parameters
# ===================================================

# -------------- Gravity Turn --------------
ALT_INITIAL_KICK = 150          # altitude to start gravity turn; [m]
DURATION_INITIAL_KICK = 30.     # duration of gravity turn; [s]


# -------------- Time Step --------------
TIME_STEP = 0.01                 # step size for integration; [s]


# -------------- Desired Orbit --------------
# ALT_DESIRED = 35786e3          # altitude of desired orbit; [m]
ALT_DESIRED = 950e3             # altitude of desired orbit; [m]
INC_DESIRED = np.deg2rad(0.)     # inclination of desired orbit; [rad]


# -------------- Launch Site --------------
LAUNCH_LAT = np.deg2rad(0.)     # latitude of launch site; [rad]
LAUNCH_LON = np.deg2rad(0.)     # longitude of launch site; [rad]


# -------------- Optimization --------------
ALPHA_LOWEST = -np.deg2rad(15.)
ALPHA_HIGHEST = -np.deg2rad(0.5)
ALPHA_INITIAL_GUESS = - np.deg2rad(10.5)


# ===================================================
# Single Run specific parameters
# ===================================================
# SS_THROTTLE = 0.35028862953186024     # Second Stage throttle
SS_THROTTLE = 1.0     # Second Stage throttle 
INITIAL_KICK_ANGLE = - np.deg2rad(2.)   # Initial kick angle [rad]


# ===================================================
# Single Run Full specific parameters
# ===================================================
DURATION_AFTER_SIMULATION = 5000.  # duration of simulation after reaching desired orbit; [s]



# ===================================================
# Direct Orbit Injection No Coast specific parameters
# ===================================================


# ===================================================
# Orbit Injection With Coasting specific parameters
# ===================================================
ALT_NO_ATMOSPHERE = 65e3        # altitude where atmosphere can be neglected; [m]


# MAX_ACCEPTED_DELTA_V = 100000.     # maximum accepted delta v; [m/s]
MAX_ACCEPTED_DELTA_V = 300.     # maximum accepted delta v; [m/s]




# ===================================================
# FOR DEBUGGING
# ===================================================

INTERRUPTS_PRINT = False
EVENTS_PRINT = False