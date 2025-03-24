""" ===============================================
                  INITIALIZATION
=============================================== """

import numpy as np

# -------------- Select Rocket ----------------
LV = 'MK1'                      # Launch Vehicle selected

# -------------- Select Simulation Type ----------------

SYM_TYPE = 4

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
ALT_INITIAL_KICK = 1000         # altitude to start gravity turn; [m]
DURATION_INITIAL_KICK = 20.     # duration of gravity turn; [s]


# -------------- Time Step --------------
TIME_STEP = 0.001                 # step size for integration; [s]


# -------------- Desired Orbit --------------
# ALT_DESIRED = 35786e3             # altitude of desired orbit; [m]
ALT_DESIRED = 300e3             # altitude of desired orbit; [m]
INC_DESIRED = np.deg2rad(0.)    # inclination of desired orbit; [rad]


# -------------- Launch Site --------------
LAUNCH_LAT = np.deg2rad(0.)     # latitude of launch site; [rad]
LAUNCH_LON = np.deg2rad(0.)     # longitude of launch site; [rad]


# -------------- Optimization --------------
ALPHA_LOWEST = -np.deg2rad(90)
ALPHA_HIGHEST = -np.deg2rad(0.1)


# ===================================================
# Single Run specific parameters
# ===================================================
# SS_THROTTLE = 0.35028862953186024     # Second Stage throttle
SS_THROTTLE = 1.0     # Second Stage throttle 
INITIAL_KICK_ANGLE = np.deg2rad(-15)   # Initial kick angle [rad]


# ===================================================
# Single Run Full specific parameters
# ===================================================
DURATION_AFTER_SIMULATION = 1000.  # duration of simulation after reaching desired orbit; [s]



# ===================================================
# Direct Orbit Injection No Coast specific parameters
# ===================================================


# ===================================================
# Orbit Injection With Coasting specific parameters
# ===================================================
ALT_NO_ATMOSPHERE = 65e3        # altitude where atmosphere can be neglected; [m]



