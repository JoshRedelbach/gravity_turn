""" ===============================================
                SOYUZ LAUNCH VEHICLE
=============================================== """

# -------------- Payload Mass --------------
M_PAYLOAD = 0e3           # payload mass; [kg]

# -------------- Event Intervals --------------
# Define time steps for events after main engine cutoff
DELAY_TIME_STAGE_SEPARATION = 3             # time when stage separation should take place after main engine cutoff 
DELAY_TIME_SECOND_ENGINE_IGNITION = 11      # time when second stage should be ignited after main engine cutoff

# -------------- Aerodynamic Properties --------------
A = 9             # cross sectional area; [m^2]
C_D = 0.3           # drag coefficient; [no unit]
C_L = 0.1           # lift coefficient; [no unit]

# =======================================================
#  FIRST STAGE
# =======================================================

# -------------- Engine Properties --------------
ISP_1 = 280             # specific impulse; [s]
F_THRUST_1 = 1800e3     # thrust of engine; [N]

# -------------- Mass Properties --------------
M_STRUCTURE_1 = 11e3   # mass structure; [kg]
M_PROP_1 = (129e3-11e3)        # mass propellant; [kg]


# =======================================================
#  SECOND STAGE
# =======================================================

# -------------- Engine Properties --------------
ISP_2 = 359              # specific impulse; [s]
F_THRUST_2 = 294e3       # thrust of engine; [N]

# -------------- Mass Properties --------------
M_STRUCTURE_2 = 2380            # mass structure; [kg]
M_PROP_2 = 25380-2380                # mass propellant; [kg]