""" ===============================================
                  INITIALIZATION
=============================================== """

import numpy as np

# -------------- Select Rocket ----------------
LV = 'MK1'                  # Launch Vehicle selected

# -------------- Select Simulation Type ----------------
SYM_TYPE = 2

"""
    1 -> Single Run
    2 -> Direct Orbit Injection No Coast
"""

# ===================================================
# General parameters
# ===================================================

# -------------- Gravity Turn --------------
alt_initial_kick = 1000                     # altitude to start gravity turn; [m]
duration_initial_kick = 20.                 # duration of gravity turn; [s]


# -------------- Time Step --------------
time_step = 0.1                   # step size for integration; [s]


# -------------- Desired Orbit --------------
alt_desired = 300e3         # altitude of desired orbit; [m]

# ===================================================
# Single Run specific parameters
# ===================================================





