""" ===============================================
                MK1 LAUNCH VEHICLE
=============================================== """

import params.constants as c

# -------------- Payload Mass --------------
m_payload = 0e3           # payload mass; [kg]

# -------------- Event Intervals --------------
# Define time steps for events after main engine cutoff
delta_time_stage_separation = 3             # time when stage separation should take place after main engine cutoff 
delta_time_second_engine_ignition = 11      # time when second stage should be ignited after main engine cutoff

# -------------- Aerodynamic Properties --------------
A = 10.             # cross sectional area; [m^2]
c_D = 0.3           # drag coefficient; [no unit]
c_L = 0.1           # lift coefficient; [no unit]

# =======================================================
#  FIRST STAGE
# =======================================================

# -------------- Engine Properties --------------
Isp_1 = 300             # specific impulse; [s]
F_thrust_1 = 4000e3     # thrust of engine; [N]

# -------------- Mass Properties --------------
m_structure_1 = 22222   # mass structure; [kg]
m_prop_1 = 200e3        # mass propellant; [kg]


# =======================================================
#  SECOND STAGE
# =======================================================

# -------------- Engine Properties --------------
Isp_2 = 420              # specific impulse; [s]
F_thrust_2 = 900e3       # thrust of engine; [N]

# -------------- Mass Properties --------------
m_structure_2 = 3333            # mass structure; [kg]
m_prop_2 = 30e3 + 12e3          # mass propellant; [kg]