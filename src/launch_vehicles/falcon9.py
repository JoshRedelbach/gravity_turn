""" ===============================================
                Falcon 9 LAUNCH VEHICLE
=============================================== """

# -------------- Payload Mass --------------
M_PAYLOAD = 0e3           # payload mass; [kg]

# -------------- Event Intervals --------------
# Define time steps for events after main engine cutoff
DELAY_TIME_STAGE_SEPARATION = 3             # time when stage separation should take place after main engine cutoff 
DELAY_TIME_SECOND_ENGINE_IGNITION = 11      # time when second stage should be ignited after main engine cutoff

# -------------- Aerodynamic Properties --------------
A = 10.52             # cross sectional area; [m^2]
C_D = 0.3           # drag coefficient; [no unit]
C_L = 0.1           # lift coefficient; [no unit]

# =======================================================
#  FIRST STAGE
# =======================================================

# -------------- Engine Properties --------------
ISP_1 = 283             # specific impulse; [s]
F_THRUST_1 = 7600e3     # thrust of engine; [N]

# -------------- Mass Properties --------------
M_STRUCTURE_1 = 25.6e3   # mass structure; [kg]
M_PROP_1 = 395.7e3        # mass propellant; [kg]


# =======================================================
#  SECOND STAGE
# =======================================================

# -------------- Engine Properties --------------
ISP_2 = 348              # specific impulse; [s]
F_THRUST_2 = 934e3       # thrust of engine; [N]

# -------------- Mass Properties --------------
M_STRUCTURE_2 = 3900            # mass structure; [kg]
M_PROP_2 = 92670                # mass propellant; [kg]