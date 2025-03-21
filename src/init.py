""" ===============================================
                  INITIALIZATION
=============================================== """

import numpy as np

# -------------- Select Rocket ----------------
LV = 'MK1'                      # Launch Vehicle selected

# -------------- Select Simulation Type ----------------
SYM_TYPE = 1
"""
    1 -> Single Run
    2 -> Direct Orbit Injection No Coast
"""

# ===================================================
# General parameters
# ===================================================

# -------------- Gravity Turn --------------
alt_initial_kick = 1000         # altitude to start gravity turn; [m]
duration_initial_kick = 20.     # duration of gravity turn; [s]


# -------------- Time Step --------------
time_step = 0.0001                 # step size for integration; [s]


# -------------- Desired Orbit --------------
alt_desired = 300e3             # altitude of desired orbit; [m]
inc_desired = np.deg2rad(0.)    # inclination of desired orbit; [rad]


# -------------- Launch Site --------------
launch_lat = np.deg2rad(0.)     # latitude of launch site; [rad]
launch_lon = np.deg2rad(0.)     # longitude of launch site; [rad]



# ===================================================
# Single Run specific parameters
# ===================================================

SS_throttle = 0.35029025077819814       # Second Stage throttle 

initial_kick_angle = np.deg2rad(-12.519912475347514) # Initial kick angle [rad]


# ===================================================
# Direct Orbit Injection No Coast specific parameters
# ===================================================




