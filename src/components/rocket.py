""" ===============================================
    Functions to simulate rocket's behaviour
=============================================== """

import components.environment as env
import components.constants as c
import init
import plotting.plot_results as plot_results
import numpy as np
from scipy.integrate import solve_ivp
from components.rocket_selector import select_rocket

# Selected Launch Vehicle
par_roc = select_rocket(init.LV)  # Replace 'MK1' with the name of your desired rocket module

#===================================================
# Interrupt functions for simulation
#===================================================

def interrupt_radius_check(t, y, SS_throttle, initial_kick_angle):
    """
    Returns zero, if the current radius exceeds the radius of the desired.
    
    Input:
        - t: current time since launch; [s]
        - y: current state vector
    """
    margin = 50e3
    r = y[1]
    if r > (init.alt_desired + c.r_earth + margin):
        #print("Interrupt Radius Check happened at time ", t)
        return 0
    return 1


def interrupt_stage_separation(t, y, SS_throttle, initial_kick_angle):
    """
    Returns zero, if the stage separation was performed successfully (everything that need to be done until the m_structure_1 needs to be separated.)
    
    Input:
        - t: current time since launch; [s]
        - y: current state vector
    """
    global time_main_engine_cutoff, main_engine_cutoff

    if main_engine_cutoff:
        if t >= (time_main_engine_cutoff + par_roc.delta_time_stage_separation):
            #print("Interrupt Stage Separation happened at time ", t)
            return 0
    return 1


# def interrupt_orbit_reached(t, y, SS_throttle, initial_kick_angle):
#     """
#     Returns zero, if the desired orbit is reached successfully meaning that current velocity norm, radius and flight path angle fit the requirements (within a certain margin).
    
#     Input:
#         - t: current time since launch; [s]
#         - y: current state vector
#     """
#     _, r, v, gamma, _ = y
    
#     r_desired = c.r_earth + init.alt_desired
#     v_desired = np.sqrt(c.mu_earth / r_desired)

#     epsilon_r = 100                        # margin for the orbit radius check
#     epsilon_v = 1                          # margin for the orbit velocity check
#     epsilon_gamma = np.deg2rad(0.1)          # margin for the flight path angle check

#     if abs(r_desired - r) < epsilon_r and abs(v_desired - v) < epsilon_v and abs(gamma) < epsilon_gamma:
#         #print("Interrupt Desired Orbit reached at time ", t)
#         return 0
#     return 1


def interrupt_stage_2_burnt(t, y, SS_throttle, initial_kick_angle):
    """
    Returns zero, if the second stage is fully burnt.
    
    Input:
        - t: current time since launch; [s]
        - y: current state vector
    """
    m = y[4]
    if m <= (par_roc.m_payload + par_roc.m_structure_2):
        #print("Interrupt Stage 2 Burnt happened at time ", t)
        return 0
    return 1


def interrupt_ground_collision(t, y, SS_throttle, initial_kick_angle):
    """
    Returns zero, if the current radius is below radius of earth
    
    Input:
        - t: current time since launch; [s]
        - y: current state vector
    """
    r = y[1]
    if r < c.r_earth - 1e3:
        #print("Interrupt Earth Collision happened at time ", t)
        return 0
    return 1


def interrupt_velocity_exceeded(t, y, SS_throttle, initial_kick_angle):
    """
    Returns zero, if the current velocity exceeds the velocity of the desired orbit.
    
    Input:
        - t: current time since launch; [s]
        - y: current state vector
    """
    v = y[2]
    r_desired = c.r_earth + init.alt_desired
    v_desired = np.sqrt(c.mu_earth / r_desired)
    margin = 0            # [m/s]
    # if v > (v_desired + margin):
    #     #print("Interrupt Desired Velocity Exceeded happened at time ", t, "\n")
    #     return 0
    return v - v_desired

    
    
#===================================================
# Event functions
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
    
    global time_main_engine_cutoff, second_engine_ignition
    
    if t >= (par_roc.delta_time_second_engine_ignition + time_main_engine_cutoff):
        second_engine_ignition = True
        
    return



#===================================================
# Define functions
#===================================================

def cartesian_coordinates(h, s):
    """
    """
    
    theta = s / c.r_earth
    y = (h + c.r_earth) * np.cos(theta)
    x = (h + c.r_earth) * np.sin(theta)
    
    return x, y

def get_orbital_elements(r, v_inertial, gamma_inertial, mu = c.mu_earth):
    """
    Computes the orbital parameters given the input state.
    
    NOTE: Note that velocity and gamma should be the relative to the inertial (ECI) reference frame, not the ECEF frame.
    
    Input:
        - r: radial distance to Earth's center; [m]
        - v_inertial: velocity relative to ECI frame; [m/s]
        - gamma_inertial: flight path angle relative to ECI frame; [rad]
    
    Output:
        - a: semimajor axis; [m]
        - e: eccentricity;
        - r_apo: apoapsis radius; [m]
        - r_peri: apoapsis radius; [m]
    """
    
    a = (mu*r) / ((2*mu) - (r*v_inertial**2))
    e = (1 - (r*v_inertial*np.cos(gamma_inertial))**2 / (mu*a))**0.5
    r_apo = a * (1+e)
    r_peri = a * (1-e)
    
    return a, e, r_apo, r_peri
    

def thrust_Isp(SS_throttle):
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
        F_T = par_roc.F_thrust_2 * SS_throttle
        Isp = par_roc.Isp_2
    else:
        print("Warning: Both first stage and second stage engines are running at the same time.")
        
    return F_T, Isp



def pitch_programm_linear(t, initial_kick_angle):
    """
    Returns the angle of attack for the initial kick.
    Increases the angle of attack to a certain value and decreases it afterwards in a linear way.

    Input:
        - t: current time since launch; [s]
    """

    global time_kick_start, kick_performed, time_raise

    if time_kick_start == None:                                     # check if kick has started
        time_kick_start = t
        #print("\nInitial kick started at t = ", t)
        return 0.0
    
    elif t > (time_kick_start + init.duration_initial_kick):     # check if kick has ended
        kick_performed = True
        #print("\nInitial kick ended at t = ", t)
        return 0.0
    
    else:
        # check if angle should raise or decrease
        if t < (time_kick_start + time_raise):
            # define rate of angle change
            angle_rate = (t - time_kick_start) / (time_raise)
            return initial_kick_angle * angle_rate
        else:
            # define rate of angle change
            angle_rate = (t - (time_kick_start + time_raise)) / (time_raise)
            return initial_kick_angle * (1 - angle_rate)



def rocket_dynamics(t, state, SS_throttle, initial_kick_angle):
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
    global time_kick_start, kick_performed, main_engine_cutoff

    # Get state components
    s, r, v, gamma, m = state

    # Compute altitude above Earth's surface
    alt = r - c.r_earth                # altitude of the rocket; [m]

    # Check main engine state and second engine state
    event_main_engine_cutoff(t, state)
    if main_engine_cutoff:
        event_second_engine_ignition(t)
    
    # --- Get current thrust, Isp ---
    F_T, Isp = thrust_Isp(SS_throttle)

    # --- Get current angle of attack ---
    if alt > init.alt_initial_kick and (not kick_performed):
        alpha = pitch_programm_linear(t, initial_kick_angle)
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
    dmdt = - F_T / (Isp * c.g0)

    return [dsdt, drdt, dvdt, dgammadt, dmdt]



def simulate_trajectory(init_time, time_stamp, state_init, stage_1_flag, stage_2_flag, SS_throttle, initial_kick_angle):
    """
    Simulates the trajectory of the rocket until a given time stamps or until a certain interrupt function is called.

    Input:
        - time_stamp: time stamp until the simulation should be performed; [s]
        - state_init: initial state vector of the rocket
    
    Output:
        - solution array of the simulation
    """

    t_span = (init_time, init_time + time_stamp + 1)
    t_eval = np.arange(init_time, init_time + time_stamp + init.time_step, init.time_step)

    if stage_1_flag:
        interrupt_list = [interrupt_radius_check, interrupt_stage_separation, interrupt_ground_collision, interrupt_velocity_exceeded]
    elif stage_2_flag:
        interrupt_list = [interrupt_radius_check, interrupt_stage_2_burnt, interrupt_ground_collision, interrupt_velocity_exceeded]
    else:
        interrupt_list = [interrupt_ground_collision]
    
    for interrupt in interrupt_list:
        interrupt.terminal = True
        interrupt.direction = 0
        
    #print(main_engine_cutoff)
    
    return solve_ivp(rocket_dynamics, y0=state_init, t_span=t_span, t_eval=t_eval, max_step=1, events=interrupt_list, atol=1e-8, args=(SS_throttle, initial_kick_angle))


# TRY TO PUT THE RUN FUNCTIONS IN THE CORRESPONDING SIM FILES <------------------------------------------------ !!!

def run(SS_throttle, initial_kick_angle):
    
    global time_kick_start, kick_performed, time_raise, main_engine_cutoff, second_engine_ignition, stage_2_burnt, time_main_engine_cutoff
    
    #===================================================
    # Reset global variables
    #===================================================
    time_kick_start = None                          # time when the initial kick starts
    kick_performed = False                          # flag to check if the initial kick has been performed
    time_raise = init.duration_initial_kick / 2. # time to raise the angle of attack to the maximum value; [s]
    main_engine_cutoff = False                      # flag to check if the first stage engine is cutoff
    second_engine_ignition = False                  # flag to check if the second stage engine is ignited
    stage_2_burnt = False                           # flag to check if the second stage is burnt
    time_main_engine_cutoff = None                  # time when the main engine cuts off

    # ---- Debugging ---- 
    # Print desired orbit
    r_desired = c.r_earth + init.alt_desired
    v_desired = np.sqrt(c.mu_earth / r_desired)
    #print(r_desired)
    #print(v_desired)

    #===================================================
    # Simulation until stage separation
    #===================================================

    # Define initial state
    initial_mass = par_roc.m_structure_1 + par_roc.m_prop_1 + par_roc.m_structure_2 + par_roc.m_prop_2 + par_roc.m_payload
    initial_state_1 = [0., c.r_earth, 0., np.deg2rad(90.), initial_mass]

    # Define time of simulation 1
    time_1 = 500.   #<------TODO

    # Call simulation for stage 1
    sol_1 = simulate_trajectory(0, time_1, initial_state_1, True, False, SS_throttle, initial_kick_angle)

    
    #===================================================
    # Simulation after stage separation
    #===================================================
    
    # Define new initial state
    initial_state_2 = sol_1.y[:, -1]

    # Adjust mass -> perform stage separation
    initial_state_2[4] = initial_state_2[4] - par_roc.m_structure_1
    
    # Define time of simulation 2
    init_time_2 = sol_1.t[-1]
    time_2 = 500.   #<------TODO
    
    # Call simulation for stage 1
    #print("Second Simulation started!")
    sol_2 = simulate_trajectory(init_time_2, time_2, initial_state_2, False, True, SS_throttle, initial_kick_angle)
    
    data = np.concatenate((sol_1.y, sol_2.y), axis=1)
    time_steps_simulation = np.concatenate((sol_1.t, sol_2.t))
    
    return time_steps_simulation, data


def run_full(SS_throttle, initial_kick_angle):
    
    global time_kick_start, kick_performed, time_raise, main_engine_cutoff, second_engine_ignition, stage_2_burnt, time_main_engine_cutoff
    
    #===================================================
    # Reset global variables
    #===================================================
    time_kick_start = None                          # time when the initial kick starts
    kick_performed = False                          # flag to check if the initial kick has been performed
    time_raise = init.duration_initial_kick / 2. # time to raise the angle of attack to the maximum value; [s]
    main_engine_cutoff = False                      # flag to check if the first stage engine is cutoff
    second_engine_ignition = False                  # flag to check if the second stage engine is ignited
    stage_2_burnt = False                           # flag to check if the second stage is burnt
    time_main_engine_cutoff = None                  # time when the main engine cuts off

    # ---- Debugging ---- 
    # Print desired orbit
    r_desired = c.r_earth + init.alt_desired
    v_desired = np.sqrt(c.mu_earth / r_desired)
    #print(r_desired)
    #print(v_desired)

    #===================================================
    # Simulation until stage separation
    #===================================================

    # Define initial state
    initial_mass = par_roc.m_structure_1 + par_roc.m_prop_1 + par_roc.m_structure_2 + par_roc.m_prop_2 + par_roc.m_payload
    initial_state_1 = [0., c.r_earth, 0., np.deg2rad(90.), initial_mass]

    # Define time of simulation 1
    time_1 = 500.   #<------TODO

    # Call simulation for stage 1
    sol_1 = simulate_trajectory(0, time_1, initial_state_1, True, False, SS_throttle, initial_kick_angle)

    
    #===================================================
    # Simulation after stage separation
    #===================================================
    
    # Define new initial state
    initial_state_2 = sol_1.y[:, -1]

    # Adjust mass -> perform stage separation
    initial_state_2[4] = initial_state_2[4] - par_roc.m_structure_1
    
    # Define time of simulation 2
    init_time_2 = sol_1.t[-1]
    time_2 = 500.   #<------TODO
    
    # Call simulation for stage 1
    #print("Second Simulation started!")
    sol_2 = simulate_trajectory(init_time_2, time_2, initial_state_2, False, True, SS_throttle, initial_kick_angle)
    
    # Cutoff second stage engine
    main_engine_cutoff = True
    second_engine_ignition = False
    
    #===================================================
    # Simulation after second engine burnout
    #===================================================
    
    # Define new initial state
    initial_state_3 = sol_2.y[:, -1]
    
    # Define time of simulation 3
    init_time_3 = sol_2.t[-1]
    time_3 = 1200.   #<------ TODO
    
    # Call simulation
    #print("Third Simulation started!")
    sol_3 = simulate_trajectory(init_time_3, time_3, initial_state_3, False, False, SS_throttle, initial_kick_angle)
    
    data = np.concatenate((sol_1.y, sol_2.y, sol_3.y), axis=1)
    time_steps_simulation = np.concatenate((sol_1.t, sol_2.t, sol_3.t))
    
    return time_steps_simulation, data
