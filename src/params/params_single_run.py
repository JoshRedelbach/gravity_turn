""" ===============================================
            PARAMETERS FOR SINGLE RUN
=============================================== """

import params.constants as c

# -------------- Rocket Properties --------------
A = 10.             # cross sectional area; [m^2]
c_D = 0.3           # drag coefficient; [no unit]
c_L = 0.1           # lift coefficient; [no unit]


# -------------- Engine Properties --------------

# We have to consider equal mass ratios for both stages!
# That means:
# m_structure_1 / (m_structure_1 + m_prop_1)
# has to be equal to
# (m_structure_2) / (m_structure_2 + m_prop_2)
# For the following values, the mass ratios are equal to 0.1.

# 1. Stage
Isp_1 = 300             # specific impulse; [s]
F_thrust_1 = 4000e3     # thrust of engine; [N]
m_structure_1 = 22222   # mass structure; [kg]
m_prop_1 = 200e3        # mass propellant; [kg]


# 2. Stage
Isp_2 = 420             # specific impulse; [s]
F_thrust_2 = 900e3      # thrust of engine; [N]
m_structure_2 = 3333    # mass structure; [kg]
m_prop_2 = 30e3         # mass propellant; [kg]


# -------------- Gravitiy Turn --------------
init_angle_of_attack = -2.      # initial angle of attack; [deg]
t_initial_kick = 20.            # time to start gravity turn; [s]
duration_initial_kick = 15.     # duration of gravity turn; [s]


# -------------- TO BE IMPLEMENTED --------------
# -------------- Desired Orbit for Final Conditions Check --------------
alt_desired = 200e3   # altitude of desired orbit; [m]
# ...

