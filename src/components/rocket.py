""" ===============================================
    Functions to simulate rocket's behaviour
=============================================== """

import components.environment as env
import params.constants as c
import params.params_rocket as par_roc
import params.params_simulation as par_sim

import numpy as np
from scipy.integrate import solve_ivp


#===================================================
# Define global variables
#===================================================
time_kick_start = None                          # time when the initial kick starts
kick_performed = False                          # flag to check if the initial kick has been performed
time_raise = par_sim.duration_initial_kick / 2. # time to raise the angle of attack to the maximum value; [s]
stage_separation = False                        # flag to check if stage separation can be performed (-> interrupt simulation 1)
separation_performed = False                    # flag to check if the stage separation has been performed
main_engine_cutoff = False                      # flag to check if the first stage engine is cutoff
second_engine_ignition = False                  # flag to check if the second stage engine is ignited 
time_main_engine_cutoff = None                  # time when the main engine cuts off


#===================================================
# Event functions to interrupt
#===================================================

def event_radius_check(t, y):
    """
    Returns zero, if the current radius exceeds the radius of the desired.
    
    Input:
        - t: current time since launch; [s]
        - y: current state vector
    """
    r = y[1]
    if r > (par_sim.alt_desired + c.r_earth):
        return 0
    else: return 1


def event_stage_separation(t, y):
    """
    Returns zero, if the stage separation was performed successfully (everything that need to be done until the m_structure_1 needs to be separated.)
    
    Input:
        - t: current time since launch; [s]
        - y: current state vector
    """
    if stage_separation: return 0
    else: return 1


def event_orbit_reached(t, y):
    """
    Returns zero, if the desired orbit is reached successfully meaning that current velocity norm, radius and flight path angle fit the requirements (within a certain margin).
    
    Input:
        - t: current time since launch; [s]
        - y: current state vector
    """
    _, r, v, gamma, _ = y
    
    r_desired = c.r_earth + par_sim.alt_desired
    v_desired = np.sqrt(c.mu_earth / r_desired)

    epsilon_r = 10              # margin for the orbit radius check
    epsilon_v = 10              # margin for the orbit velocity check
    epsilon_gamma = 10          # margin for the flight path angle check

    if abs(r_desired - r) < epsilon_r and abs(v_desired - v) < epsilon_v and abs(gamma) < epsilon_gamma:
        return 0
    else: return 1
    
    
#===================================================
# Event functions to interrupt
#===================================================

def event_main_engine_cutoff(t, y):
    """
    Checks if there is still propellant in the first stage and triggers engine cutoff event when there isn't.
    
    Input:
        - y: current state vector
    """
    
    global main_engine_cutoff, time_main_engine_cutoff
    
    if main_engine_cutoff == True:
        return
    
    first_stage_leftover_propellant = y[4] - (par_roc.m_structure_1 + par_roc.m_structure_2 + par_roc.m_prop_2 + par_roc.m_payload)
    
    if first_stage_leftover_propellant <= 0 and main_engine_cutoff == False:
        main_engine_cutoff = True
        time_main_engine_cutoff = t

    return

def event_second_engine_ignition(t):
    """
    Triggers second stage engine ignition.
    
    Input:
        - t: current time since launch; [s]
    """
    
    global delta_time_second_engine_ignition, time_main_engine_cutoff, second_engine_ignition
    
    if t >= (delta_time_second_engine_ignition + time_main_engine_cutoff):
        second_engine_ignition = True
        
    return
    
    
    



#===================================================
# Define functions
#===================================================

def thrust_Isp():
    """
    Returns the current thrust and Isp.
    
    NOTE: It doesn't account for both engines ignited.
    
    Output:
        - F_T: current thrust; [N]
        - Isp: current specific impulse; [s]
    """
    
    global main_engine_cutoff, second_engine_ignition
    
    if main_engine_cutoff == False:
        F_T = par_roc.F_thrust_1
        Isp = par_roc.Isp_1
    elif main_engine_cutoff == True and second_engine_ignition == False:
        F_T = 0
        Isp = par_roc.Isp_1
    elif main_engine_cutoff == True and second_engine_ignition == True:
        F_T = par_roc.F_thrust_2
        Isp = par_roc.Isp_2
    else:
        print("Warning: Both first stage and second stage engines are running at the same time.")
        
    return F_T, Isp



def pitch_programm_linear(t):
    """
    Returns the angle of attack for the initial kick.
    Increases the angle of attack to a certain value and decreases it afterwards in a linear way.

    Input:
        - t: current time since launch; [s]
    """

    global time_kick_start, kick_performed, time_raise

    if time_kick_start == None:                                     # check if kick has started
        time_kick_start = t
        print("\nInitial kick started at t = ", t)
        return 0.0
    
    elif t > (time_kick_start + par_sim.duration_initial_kick):     # check if kick has ended
        kick_performed = True
        print("\nInitial kick ended at t = ", t)
        return 0.0
    
    else:
        # check if angle should raise or decrease
        if t < (time_kick_start + time_raise):
            # define rate of angle change
            angle_rate = (t - time_kick_start) / (time_raise)
            return par_sim.max_angle_of_attack * angle_rate
        else:
            # define rate of angle change
            angle_rate = (t - (time_kick_start + time_raise)) / (time_raise)
            return par_sim.max_angle_of_attack * (1 - angle_rate)



def rocket_dynamics(t, state):
    """
    Simulates the dynamics of the rocket. This function will be integrated by the scipy.solve_ivp function
    
    Input:
        - t: time variable (necessary for solve_ivp function)
        - state: current state vector of the rocket

    NOTE: The state vector is defined as follows:
        - s: downtrack; [m]
        - r: radius from Earth's center; [m]
        - v: velocity norm; [m/s]
        - gamma: flight path angle; [rad]
        - m: current mass; [kg]

    Output:
        - derivatives of the state vector
    """
    global time_kick_start, kick_performed

    # Get state components
    s, r, v, gamma, m = state

    # Compute altitude above Earth's surface
    alt = r - c.r_earth                # altitude of the rocket; [m]

    # Check main engine state and second engine state
    event_main_engine_cutoff(t, state)
    event_second_engine_ignition(t)
    
    # --- Get current thrust, Isp ---
    F_T, Isp = thrust_Isp()

    # --- Get current angle of attack ---
    if alt > par_sim.alt_initial_kick and (not kick_performed):
        alpha = pitch_programm_linear(t)
    else:
        alpha = 0.

    # --- Determine current accelerations and forces ---
    a_grav = env.accel_grav(r)                                  # gravity at the altitude of the rocket in negative radial direction
    F_D = env.drag_force(v, alt, par_roc.c_D, par_roc.A)        # drag force norm acting in negative velocity direction
    
    # NOTE: not sure, if we have to include lift forces or not
    # F_L = env.lift_force(v, alt, par_roc.c_L, par_roc.A)      # lift force norm acting in vertical to velocity direction
    F_L = 0.0

    # --- Get trigonometric operations of gamma and alpha ---
    c_gamma = np.cos(gamma)
    s_gamma = np.sin(gamma)
    c_alpha = np.cos(alpha)
    s_alpha = np.sin(alpha)

    # --- Compute the derivatives ---
    dsdt = (c.r_earth / r) * v * c_gamma

    # Reference: [1] Equations 6.3.8
    drdt = v * np.sin(gamma)

    # Reference: [1] Equation 6.3.6
    dvdt = (F_T / m) * c_alpha - (F_D / m) - a_grav * s_gamma

    # Catch the case of zero velocity to avoid division by zero
    epsilon = 1e-6
    if v < epsilon:
        dgammadt = 0.
    else:
        # Reference: [1] Equation 6.3.7
        dgammadt = (1./v) * ( (F_T / m) * s_alpha + F_L / m  - (a_grav - (v**2 / r)) * c_gamma )
    
    # Derivative of mass
    if t > par_roc.t_burn_1 and (m > (par_roc.m_structure_2 + par_roc.m_prop_2)):
        # If first stage is burnt out, and the mass still includes the mass of the first stage structure, subtract it ONCE!
        dmdt = - (m - par_roc.m_structure_1) # * ( (par_roc.t_burn_1 + par_roc.t_burn_2) / par_sim.number_of_points )
        print("\nFirst stage burnt out at t = ", t)
        print("\nStage separation performed. Current mass: ", m)
    else: dmdt = - F_T / (Isp * c.g0)

    return [dsdt, drdt, dvdt, dgammadt, dmdt]



def simulate_trajectory(time_stamp, state_init):
    """
    Simulates the trajectory of the rocket until a given time stamps.

    Input:
        - time_stamp: time stamp until the simulation should be performed; [s]
        - state_init: initial state vector of the rocket
    
    Output:
        - solution array of the simulation
    """

    t_span = (0, time_stamp)
    t_eval = np.arange(0.0, time_stamp + par_sim.time_step, par_sim.time_step)
    
    return solve_ivp(rocket_dynamics, y0=state_init, t_span=t_span, t_eval=t_eval, max_step=1)