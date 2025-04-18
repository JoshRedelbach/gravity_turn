""" ===============================================
    Functions to simulate rocket's behaviour
=============================================== """

import components.environment as env
import components.constants as c
import init
import numpy as np
from scipy.integrate import solve_ivp
from components.rocket_selector import select_rocket
import components.solvers as solvers
import simulation.coasting_single_burn as coasting_single_burn
import simulation.coasting_double_burn as coasting_double_burn

# Selected Launch Vehicle
par_roc = select_rocket(init.LV)

#===================================================
# Interrupt functions for simulation
#===================================================

def interrupt_radius_check(t, y, ss_throttle, initial_kick_angle):
    """
    Returns zero, if the current radius exceeds the radius of the desired.
    
    Input:
        - t: current time since launch; [s]
        - y: current state vector
        - ss_throttle: throttle of the second stage
        - initial_kick_angle: angle of attack for the initial kick
    """
    margin = 50e3
    r = y[1]
    if r > (init.ALT_DESIRED + c.R_EARTH + margin):
        if init.INTERRUPTS_PRINT: print("Interrupt Radius Check happened at time ", t)
        return 0
    return 1



def interrupt_stage_separation(t, y, ss_throttle, initial_kick_angle):
    """
    Returns zero, if the stage separation was performed successfully (everything that need to be done until the M_STRUCTURE_1 needs to be separated.)
    
    Input:
        - t: current time since launch; [s]
        - y: current state vector
        - ss_throttle: throttle of the second stage
        - initial_kick_angle: angle of attack for the initial kick
    """
    global time_main_engine_cutoff, main_engine_cutoff

    if main_engine_cutoff:
        if t >= (time_main_engine_cutoff + par_roc.DELAY_TIME_STAGE_SEPARATION):
            if init.INTERRUPTS_PRINT: print("Interrupt Stage Separation happened at time ", t)
            return 0
    return 1


# def interrupt_orbit_reached(t, y, ss_throttle, initial_kick_angle):
#     """
#     Returns zero, if the desired orbit is reached successfully meaning that current velocity norm, radius and flight path angle fit the requirements (within a certain margin).
    
#     Input:
#         - t: current time since launch; [s]
#         - y: current state vector
#         - ss_throttle: throttle of the second stage
#         - initial_kick_angle: angle of attack for the initial kick
#     """
#     _, r, v, gamma, _ = y
    
#     r_desired = c.R_EARTH + init.ALT_DESIRED
#     v_desired = np.sqrt(c.MU_EARTH / r_desired)

#     epsilon_r = 100                        # margin for the orbit radius check
#     epsilon_v = 1                          # margin for the orbit velocity check
#     epsilon_gamma = np.deg2rad(0.1)          # margin for the flight path angle check

#     if abs(r_desired - r) < epsilon_r and abs(v_desired - v) < epsilon_v and abs(gamma) < epsilon_gamma:
#         #print("Interrupt Desired Orbit reached at time ", t)
#         return 0
#     return 1


def interrupt_stage_2_burnt(t, y, ss_throttle, initial_kick_angle):
    """
    Returns zero, if the second stage is fully burnt.
    
    Input:
        - t: current time since launch; [s]
        - y: current state vector
        - ss_throttle: throttle of the second stage
        - initial_kick_angle: angle of attack for the initial kick
    """
    m = y[4]
    if m <= (par_roc.M_PAYLOAD + par_roc.M_STRUCTURE_2):
        if init.INTERRUPTS_PRINT: print("Interrupt Stage 2 Burnt happened at time ", t)
        return 0
    return 1


def interrupt_ground_collision(t, y, ss_throttle, initial_kick_angle):
    """
    Returns zero, if the current radius is below radius of earth
    
    Input:
        - t: current time since launch; [s]
        - y: current state vector
        - ss_throttle: throttle of the second stage
        - initial_kick_angle: angle of attack for the initial kick
    """
    r = y[1]
    if r < c.R_EARTH - 1e3:
        if init.INTERRUPTS_PRINT: print("Interrupt Earth Collision happened at time ", t)
        return 0
    return 1


def interrupt_velocity_exceeded(t, y, ss_throttle, initial_kick_angle):
    """
    Returns zero, if the current velocity exceeds the velocity of the desired orbit.
    
    Input:
        - t: current time since launch; [s]
        - y: current state vector
        - ss_throttle: throttle of the second stage
        - initial_kick_angle: angle of attack for the initial kick
    """
    v = y[2]
    r_desired = c.R_EARTH + init.ALT_DESIRED
    v_desired = np.sqrt(c.MU_EARTH / r_desired)
    margin = 0            # [m/s]
    # if v > (v_desired + margin):
    #     #print("Interrupt Desired Velocity Exceeded happened at time ", t, "\n")
    #     return 0
    return v - v_desired


def interrupt_single_burn_traj(t, y, ss_throttle, initial_kick_angle):
    """
    Performed as soon as the rocket reached the altitude where the atmosphere can be neglected.
    Checks if the current apogee is within a certain margin close to the desired altitude (return 0) or not (return 1) if the rocket would stop burnig at the current time stamp.
    
    Input:
        - t: current time since launch; [s]
        - y: current state vector
        - ss_throttle: throttle of the second stage
        - initial_kick_angle: angle of attack for the initial kick
    """
    r = y[1]
    v = y[2]
    gamma = y[3]
    alt = r - c.R_EARTH

    if alt < init.ALT_NO_ATMOSPHERE:
        return 1
    else:
        # Compute current orbital elements
        a, e, r_apo, r_peri, _ = get_orbital_elements(r, v, gamma)

        diff = r_apo - (init.ALT_DESIRED + c.R_EARTH)
        
        return diff
    

def interrupt_horizontal_check(t, y, ss_throttle, initial_kick_angle):
    """
    Checks if the rocket reached horizontal flight direction.
    
    Input:
        - t: current time since launch; [s]
        - y: current state vector
        - ss_throttle: throttle of the second stage
        - initial_kick_angle: angle of attack for the initial kick
    """
    gamma = y[3]
    epsilon = np.deg2rad(0.01)
    
    if gamma < epsilon:
        if init.INTERRUPTS_PRINT: print("Interrupt Horizontal Flight Direction happened at time ", t)
        return 0
    return 1

    
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
    
    first_stage_leftover_propellant = y[4] - (par_roc.M_STRUCTURE_1 + par_roc.M_STRUCTURE_2 + par_roc.M_PROP_2 + par_roc.M_PAYLOAD)
    
    if first_stage_leftover_propellant <= 0 and main_engine_cutoff == False:
        main_engine_cutoff = True
        time_main_engine_cutoff = t

        if init.EVENTS_PRINT: print("Main engine cutoff at t = ", t)

    return


def event_second_engine_ignition(t):
    """
    Triggers second stage engine ignition.
    
    Input:
        - t: current time since launch; [s]
    """
    
    global time_main_engine_cutoff, second_engine_ignition
    
    if t >= (par_roc.DELAY_TIME_SECOND_ENGINE_IGNITION + time_main_engine_cutoff):
        second_engine_ignition = True
        
        if init.EVENTS_PRINT: print("Second engine ignited at t = ", t)
        
    return

#===================================================
# Define functions
#===================================================

def cartesian_coordinates(h, s):
    """
    """
    
    theta = s / c.R_EARTH
    y = (h + c.R_EARTH) * np.cos(theta)
    x = (h + c.R_EARTH) * np.sin(theta)
    
    return x, y

def get_orbital_elements(r, v_inertial, gamma_inertial, mu = c.MU_EARTH):
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
        - orbit_period: period of the orbit; [s]
    """
    
    a = (mu*r) / ((2*mu) - (r*v_inertial**2))
    e = (1 - (r*v_inertial*np.cos(gamma_inertial))**2 / (mu*a))**0.5
    r_apo = a * (1+e)
    r_peri = a * (1-e)
    orbit_period = 2 * np.pi * (np.pow(a, 1.5)) / (np.pow(mu, 0.5))
    
    return a, e, r_apo, r_peri, orbit_period
    


def thrust_Isp(ss_throttle):
    """
    Returns the current thrust and Isp.
    
    NOTE: It doesn't account for both engines ignited.
    
    Output:
        - F_T: current thrust; [N]
        - Isp: current specific impulse; [s]
    """
    
    global main_engine_cutoff, second_engine_ignition, second_stage_cutoff
    
    if not main_engine_cutoff:
        F_T = par_roc.F_THRUST_1
        Isp = par_roc.ISP_1
    elif main_engine_cutoff and not second_engine_ignition:
        F_T = 0
        Isp = par_roc.ISP_1
    elif main_engine_cutoff and second_stage_cutoff:
        F_T = 0
        Isp = par_roc.ISP_2
    elif main_engine_cutoff and second_engine_ignition:
        F_T = par_roc.F_THRUST_2 * ss_throttle
        Isp = par_roc.ISP_2
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
        if init.EVENTS_PRINT: print("\nInitial kick started at t = ", t)
        return 0.0
    
    elif t > (time_kick_start + init.DURATION_INITIAL_KICK):     # check if kick has ended
        kick_performed = True
        if init.EVENTS_PRINT: print("\nInitial kick ended at t = ", t)
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
        


def compute_double_burn_delta_v(state):
    """
    Computes the delta-v required to reach the desired orbit with a coasting phase and Hohman transfer.

    Input:
        - state: current state vector of the rocket
    
    Output:
        - delta_v: delta-v required to reach the desired orbit; [m/s]
    """
    # Getting the state components
    r = state[1]
    v = state[2]
    gamma = state[3]
    alt = r - c.R_EARTH
    r_desired = c.R_EARTH + init.ALT_DESIRED

    # Compute current orbital elements
    a, e, r_apo, r_peri, _ = get_orbital_elements(r, v, gamma)

    if (r_apo > r_desired) or (alt < init.ALT_NO_ATMOSPHERE):
        return 99999999., 99999999., 99999999.
    
    # Compute velocity at apogee
    v_apo = np.sqrt(c.MU_EARTH * a * (1 - e**2)) / r_apo

    # Compute delta-v for Hohman transfer
    delta_v_total, delta_v1, delta_v2 = solvers.hohman_transfer(v_apo, r_apo, r_desired)

    return delta_v_total, delta_v1, delta_v2


def check_dual_burn_prop_used(delta_v_total, delta_v1, delta_v2, state, t):
    """
    Checks if the current prop used is smaller than current smallest prop for a specific kick angle. If so, update the current smallest prop, the corresponding altitude, delta-v and time stamp.
    
    Input:
        - delta_v: computed delta-v
        - state: current state vector of the rocket
        - t: current time since launch; [s]
    """
    global double_burn_smallest_prop

    # Calculate current total propellant used of second stage
    m_propellant_left = state[4] - (par_roc.M_STRUCTURE_2 + par_roc.M_PAYLOAD)
    m_propellant_used = par_roc.M_PROP_2 - m_propellant_left
    m_propellant_required = state[4] * (1 - np.exp(-delta_v_total / (c.G0 * par_roc.ISP_2)))

    # Check if the propellant required is less than the propellant left
    if m_propellant_required < m_propellant_left:
        m_propellant_total_used_2nd_stage = m_propellant_used + m_propellant_required
    else:
        m_propellant_total_used_2nd_stage = 999999999.

    # Update best run for current kick angle if necessary
    if m_propellant_total_used_2nd_stage < double_burn_smallest_prop:
        double_burn_smallest_prop = m_propellant_total_used_2nd_stage

    # Save results of best run of whole optimization if necessary
    if m_propellant_total_used_2nd_stage < coasting_double_burn.double_burn_smallest_prop_global:
        coasting_double_burn.double_burn_altitude_global = state[1] - c.R_EARTH
        coasting_double_burn.double_burn_time_global = t
        coasting_double_burn.double_burn_smallest_prop_global = m_propellant_total_used_2nd_stage
        coasting_double_burn.double_burn_delta_v_1_global = delta_v1
        coasting_double_burn.double_burn_delta_v_2_global = delta_v2


#===================================================
# Define Dynamics
#===================================================


def rocket_dynamics(t, state, ss_throttle, initial_kick_angle):
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
        - lat: latitude (ECI); [rad]
        - lon: longitude (ECI); [rad]
        - ceta: heading angle (ECI); [rad]

    Output:
        - derivatives of the state vector
    """
    global time_kick_start, kick_performed, main_engine_cutoff, flag_falling_single_burn

    # Get state components
    s, r, v, gamma, m, lat, lon, ceta = state

    # Compute altitude above Earth's surface
    alt = r - c.R_EARTH                # altitude of the rocket; [m]

    # Check main engine state and second engine state
    event_main_engine_cutoff(t, state)
    if main_engine_cutoff:
        event_second_engine_ignition(t)
    
    # --- Get current thrust, Isp ---
    F_T, Isp = thrust_Isp(ss_throttle)

    # --- Get current angle of attack ---
    if t >= init.TIME_TO_START_KICK and (not kick_performed):
        alpha = pitch_programm_linear(t, initial_kick_angle)
    else:
        alpha = 0.

    # --- Determine current accelerations and forces ---
    a_grav = env.accel_grav(r)                                  # gravity at the altitude of the rocket in negative radial direction
    F_D = env.drag_force(v, alt, par_roc.C_D, par_roc.A)        # drag force norm acting in negative velocity direction
    
    # NOTE: not sure, if we have to include lift forces or not
    # F_L = env.lift_force(v, alt, par_roc.C_L, par_roc.A)      # lift force norm acting in vertical to velocity direction
    F_L = 0.0

    state_differentiated = diff_eom_base(s, r, v, gamma, m, F_L, F_D, F_T, a_grav, alpha, Isp)
    # state_differentiated = diff_eom_advanced(s, r, v, gamma, m, lat, lon, ceta, F_L, F_D, F_T, a_grav, alpha, Isp)

    if time_kick_start == None:
        state_differentiated[3] = 0.0

    if state_differentiated[2] < 0:
        flag_falling_single_burn = True

    # DEBUGGING
    # if t < 0.001:
    #     print("\n")
    #     print("Downtrack: ", state[0])
    #     print("Radius: ", state[1])
    #     print("Velocity: ", state[2])
    #     print("Flight path angle: ", state[3])
    #     print("Mass: ", state[4])
    #     print("Downtrack Dot: ", state_differentiated[0])
    #     print("Radius Dot: ", state_differentiated[1])
    #     print("Velocity Dot: ", state_differentiated[2])
    #     print("Flight path angle Dot: ", state_differentiated[3])
    #     print("Mass Dot: ", state_differentiated[4])
    #     print("a_grav", a_grav)
    #     print("F_D", F_D)
    #     print("F_T", F_T)
    #     print("F_L", F_L)
    #     print("alpha", alpha)
    #     print("\n")

    # ---- Check for coasting phase double burn ----
    if init.SYM_TYPE == 5 and not(coasting_double_burn.DOUBLE_BURN_FULL_SIMULATION):
        if main_engine_cutoff and second_engine_ignition:
            if state_differentiated[2] > 0:
                # Calculate delta-v for double burn
                delta_v_total, delta_v1, delta_v2 = compute_double_burn_delta_v(state)

                # if (delta_v1 < init.MAX_ACCEPTED_DELTA_V) and (delta_v2 < init.MAX_ACCEPTED_DELTA_V):
                #     # Check if the current prop used is smaller than current smallest prop for a specific kick angle
                #     check_dual_burn_prop_used(delta_v_total, delta_v1, delta_v2, state, t)

                # check if burn time of those maneuver would be realistic
                # get burn time of delta_v1 and delta_v2
                m1 = m
                m_propellant_required_1 = m1 * (1 - np.exp(-delta_v1 / (c.G0 * par_roc.ISP_2)))
                m2 = m - m_propellant_required_1
                burn_time_delta_v_1 = solvers.calculate_burn_time(m1, delta_v1)
                burn_time_delta_v_2 = solvers.calculate_burn_time(m2, delta_v2)

                if (burn_time_delta_v_1 < init.MAX_ACCEPTED_BURN_TIME) and (burn_time_delta_v_2 < init.MAX_ACCEPTED_BURN_TIME):
                    # Check if the current prop used is smaller than current smallest prop for a specific kick angle
                    check_dual_burn_prop_used(delta_v_total, delta_v1, delta_v2, state, t)


    return state_differentiated



def diff_eom_base(s, r, v, gamma, m, F_L, F_D, F_T, a_grav, alpha, Isp):
    """
    Differential equations of motion for the rocket WITHOUT earth rotation.
    
    Input:
        - s: downtrack; [m]
        - r: radius from Earth's center; [m]
        - v: velocity norm; [m/s]
        - gamma: flight path angle; [rad]
        - m: current mass; [kg]
        - F_L: lift force; [N]
        - F_D: drag force; [N]
        - F_T: thrust force; [N]
        - a_grav: gravity acceleration; [m/s^2]
        - alpha: angle of attack; [rad]
        - Isp: specific impulse; [s]

    Output:
        - derivatives of the state vector (for all variables)
    """
    # --- Get trigonometric operations of gamma and alpha ---
    c_gamma = np.cos(gamma)
    s_gamma = np.sin(gamma)
    c_alpha = np.cos(alpha)
    s_alpha = np.sin(alpha)

    # --- Compute the derivatives ---
    dsdt = (c.R_EARTH / r) * v * c_gamma

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
    dmdt = - F_T / (Isp * c.G0)

    return [dsdt, drdt, dvdt, dgammadt, dmdt, 0, 0, 0]


def diff_eom_advanced(s, r, v, gamma, m, lat, lon, ceta, F_L, F_D, F_T, a_grav, alpha, Isp):
    """ 
    NOTE: NOT TESTED YET!
    Be careful with the advanced equations of motion. They are not fully implemented yet!
    Division by zero could occur!
    """
    """
    Differential equations of motion for the rocket WITH earth rotation. 
    
    Reference: "Particle Swarm Optimization of Ascent Trajectories of Multistage Launch Vehicles (Mauro Pontani, 2013)", page 8
    
    Input:
        - s: downtrack; [m]
        - r: radius from Earth's center; [m]
        - v: velocity norm; [m/s]
        - gamma: flight path angle; [rad]
        - m: current mass; [kg]
        - lat: latitude (ECI); [rad]
        - lon: longitude (ECI); [rad]
        - ceta: heading angle (ECI); [rad]
        - F_L: lift force; [N]
        - F_D: drag force; [N]
        - F_T: thrust force; [N]
        - a_grav: gravity acceleration; [m/s^2]
        - alpha: angle of attack; [rad]
        - Isp: specific impulse; [s]

    Output:
        - derivatives of the state vector (for all variables)
    """
    # --- Get trigonometric operations of gamma and alpha ---
    c_gamma = np.cos(gamma)
    s_gamma = np.sin(gamma)
    c_alpha = np.cos(alpha)
    s_alpha = np.sin(alpha)

    # --- Compute the derivatives ---
    dsdt = (c.R_EARTH / r) * v * c_gamma

    # Reference: [1] Equations 6.3.8
    drdt = v * np.sin(gamma)

    # dvdt
    dvdt = (F_T / m) * c_alpha - (F_D / m) - a_grav * s_gamma + c.OMEGA_EARTH**2 * r * np.cos(lat) * (np.cos(lat) * np.sin(gamma) - np.sin(lat) * np.cos(gamma) * np.sin(ceta))

    # Catch the case of zero velocity to avoid division by zero
    epsilon = 1e-6
    if v < epsilon:
        dgammadt = 0.
    else:

        # 1st term
        term_1 = (F_T / m) * s_alpha

        # 2nd term
        term_2 = (v**2 * r - c.MU_EARTH) / (r**2 * v) * c_gamma

        # 3rd term
        term_3 = 2 * c.OMEGA_EARTH * np.cos(lat) * np.cos(ceta)

        # 4th term
        term_4 = ((c.OMEGA_EARTH**2 * r) / v) * np.cos(lat) * (np.cos(lat) * np.cos(gamma) + np.sin(lat) * np.sin(gamma) * np.sin(ceta))

        dgammadt = term_1 + term_2 + term_3 + term_4
    
    # Derivative of mass
    dmdt = - F_T / (Isp * c.G0)

    # Derivative of latitude
    dlatdt = (v * np.cos(gamma) * np.sin(ceta)) / r

    # Derivative of longitude
    dlondt = (v * np.cos(gamma) * np.cos(ceta)) / r

    # Derivative of heading angle ceta
    # NOTE: how do we handle the division by zero because of division by cos(gamma) and the division by zero because of v = 0?
    dcetadt = - (v /r) * np.tan(lat) * np.cos(gamma) * np.cos(ceta) + 2 * c.OMEGA_EARTH * np.cos(lat) * np.tan(gamma) * np.sin(ceta) - ((c.OMEGA_EARTH**2 * r) / (v * np.cos(gamma))) * np.sin(lat) * np.cos(lat) * np.cos(ceta) - 2 * c.OMEGA_EARTH * np.sin(lat) + 0.

    return [dsdt, drdt, dvdt, dgammadt, dmdt, dlatdt, dlondt, dcetadt]




#===================================================
# Define Simulations
#===================================================

def simulate_trajectory(init_time, time_stamp, state_init, stage_1_flag, stage_2_flag, ss_throttle, initial_kick_angle):
    """
    Simulates the trajectory of the rocket until a given time stamps or until a certain interrupt function is called.

    Input:
        - time_stamp: time stamp until the simulation should be performed; [s]
        - state_init: initial state vector of the rocket
    
    Output:
        - solution array of the simulation
    """

    t_span = (init_time, init_time + time_stamp + 1)
    t_eval = np.arange(init_time, init_time + time_stamp + init.TIME_STEP, init.TIME_STEP)

    if stage_1_flag:
        interrupt_list = [interrupt_stage_separation, interrupt_ground_collision, interrupt_velocity_exceeded]
    
    elif stage_2_flag:
        if (init.SYM_TYPE == 4):
            interrupt_list = [interrupt_radius_check, interrupt_stage_2_burnt, interrupt_ground_collision, interrupt_single_burn_traj, interrupt_horizontal_check]
        elif init.SYM_TYPE == 5 and not(coasting_double_burn.DOUBLE_BURN_FULL_SIMULATION):
            interrupt_list = [interrupt_radius_check, interrupt_stage_2_burnt, interrupt_ground_collision, interrupt_single_burn_traj, interrupt_horizontal_check]
        elif init.SYM_TYPE == 5 and (coasting_double_burn.DOUBLE_BURN_FULL_SIMULATION):
            # interrupt_list = [interrupt_radius_check, interrupt_stage_2_burnt, interrupt_ground_collision, interrupt_velocity_exceeded]
            interrupt_list = []
        else:
            interrupt_list = [interrupt_radius_check, interrupt_stage_2_burnt, interrupt_ground_collision, interrupt_velocity_exceeded]

    else:
        interrupt_list = [interrupt_ground_collision]
    
    for interrupt in interrupt_list:
        interrupt.terminal = True
        interrupt.direction = 0
    
    return solve_ivp(rocket_dynamics, y0=state_init, t_span=t_span, t_eval=t_eval, max_step=1, events=interrupt_list, atol=1e-8, args=(ss_throttle, initial_kick_angle))



# TRY TO PUT THE RUN FUNCTIONS IN THE CORRESPONDING SIM FILES <------------------------------------------------ !!!

def run(ss_throttle, initial_kick_angle):
    
    global time_kick_start, kick_performed, time_raise, main_engine_cutoff, second_engine_ignition, stage_2_burnt, time_main_engine_cutoff, second_stage_cutoff, double_burn_smallest_prop, flag_falling_single_burn
    
    #===================================================
    # Reset global variables
    #===================================================
    time_kick_start = None                          # time when the initial kick starts
    kick_performed = False                          # flag to check if the initial kick has been performed
    time_raise = init.DURATION_INITIAL_KICK / 2. # time to raise the angle of attack to the maximum value; [s]
    main_engine_cutoff = False                      # flag to check if the first stage engine is cutoff
    second_engine_ignition = False                  # flag to check if the second stage engine is ignited
    stage_2_burnt = False                           # flag to check if the second stage is burnt
    time_main_engine_cutoff = None                  # time when the main engine cuts off
    second_stage_cutoff = False                     # flag to check if the second stage is cutoff
    double_burn_smallest_prop = 9e12                # current smallest prop used for double burn for current tested kick angle
    flag_falling_single_burn = False                # flag to check if the rocket is falling back to Earth; only used for single burn optimization

    # ---- Debugging ---- 
    # Print desired orbit
    # r_desired = c.R_EARTH + init.ALT_DESIRED
    # v_desired = np.sqrt(c.MU_EARTH / r_desired)
    #print(r_desired)
    #print(v_desired)

    #===================================================
    # Simulation until stage separation
    #===================================================

    # Define initial state
    initial_mass = par_roc.M_STRUCTURE_1 + par_roc.M_PROP_1 + par_roc.M_STRUCTURE_2 + par_roc.M_PROP_2 + par_roc.M_PAYLOAD
    initial_state_1 = [0., c.R_EARTH, 0., np.deg2rad(90.), initial_mass, 0, 0, 0]

    # Define time of simulation 1
    time_1 = 500.   #<------TODO

    # Call simulation for stage 1
    sol_1 = simulate_trajectory(0, time_1, initial_state_1, True, False, ss_throttle, initial_kick_angle)

    
    #===================================================
    # Simulation after stage separation
    #===================================================

    if (init.SYM_TYPE == 5) and (coasting_double_burn.DOUBLE_BURN_FULL_SIMULATION):
        return simulate_full_double_burn_trajectory(sol_1, initial_kick_angle)
    
    # Define new initial state
    initial_state_2 = sol_1.y[:, -1]

    # Adjust mass -> perform stage separation
    initial_state_2[4] = initial_state_2[4] - par_roc.M_STRUCTURE_1
    
    # Define time of simulation 2
    init_time_2 = sol_1.t[-1]
    time_2 = 4000.   #<------TODO
    
    # Call simulation for stage 1
    # print("Second Simulation started!")
    sol_2 = simulate_trajectory(init_time_2, time_2, initial_state_2, False, True, ss_throttle, initial_kick_angle)

    if init.SYM_TYPE != 2:
        data = np.concatenate((sol_1.y, sol_2.y), axis=1)
        time_steps_simulation = np.concatenate((sol_1.t, sol_2.t))

        if init.SYM_TYPE == 4:
            
            r_stop = sol_2.y[1, -1]
            v_stop = sol_2.y[2, -1]
            gamma_stop = sol_2.y[3, -1]
        
            # Calculate altitude to stop burning
            alt_stop = r_stop - c.R_EARTH
            
            # Calculate orbital elements at stop
            a_stop, e_stop, r_apo_stop, r_peri_stop, orbit_period_stop = get_orbital_elements(r_stop, v_stop, gamma_stop)

            epsilon = (c.R_EARTH + init.ALT_DESIRED) * 0.002   # meters
            diff = abs(r_apo_stop - (c.R_EARTH + init.ALT_DESIRED))

            if diff < epsilon:
                
                # ----- Calculate delta v -----
                r_desired = c.R_EARTH + init.ALT_DESIRED
                v_desired = np.sqrt(c.MU_EARTH / r_desired)
                
                # Get velocity at apogee
                v_apo = np.sqrt(c.MU_EARTH * a_stop * (1 - e_stop**2)) / r_apo_stop
                # v_apo = np.sqrt(c.MU_EARTH * ((2. / r_apo_stop) - (1. / a_stop)))
                # ---> Josh checked, same result so both formulas are correct

                delta_v = np.abs(v_apo - v_desired)

                # ----- Calculate total propellant required -----
                m_propellant_left = sol_2.y[4, -1] - (par_roc.M_STRUCTURE_2 + par_roc.M_PAYLOAD)
                m_propellant_used = par_roc.M_PROP_2 - m_propellant_left
                m_propellant_required = sol_2.y[4, -1] * (1 - np.exp(-delta_v / (c.G0 * par_roc.ISP_2)))

                # Check if the propellant required is less than the propellant left
                if m_propellant_required < m_propellant_left:
                    m_propellant_total_used_2nd_stage = m_propellant_used + m_propellant_required
                else:
                    m_propellant_total_used_2nd_stage = 999999999.

                # Calculate burn time for that delta_v
                m1 = sol_2.y[4, -1]
                burn_time_delta_v = solvers.calculate_burn_time(m1, delta_v)
                if burn_time_delta_v > init.MAX_ACCEPTED_BURN_TIME:
                    m_propellant_total_used_2nd_stage = 999999999.

                # if delta_v > init.MAX_ACCEPTED_DELTA_V:
                #     m_propellant_total_used_2nd_stage = 999999999.

                # # Print masses for debugging
                # print("\n\n")
                # print("Propellant left: \t\t\t", m_propellant_left)
                # print("Propellant used: \t\t\t", m_propellant_used)
                # print("Propellant required by circularization: ", m_propellant_required)
                # print("Total propellant used: \t\t\t", m_propellant_total_used_2nd_stage)
                # print("\n")

                if not(coasting_single_burn.SINGLE_BURN_FULL_SIMULATION):
                    return time_steps_simulation, data, alt_stop, delta_v, m_propellant_total_used_2nd_stage
                else:
                    coasting_single_burn.TIME_TO_STOP_BURNING_SINGLE_BURN_FINAL = sol_2.t[-1]
                    par_roc.C_D = 0.
                    # Print result of masses
                    print("\t* Optimal altitude to stop burning: \t\t", alt_stop/1000, "km")
                    print("\t* Optimal time to stop burning: \t\t", coasting_single_burn.TIME_TO_STOP_BURNING_SINGLE_BURN_FINAL, "s")
                    print("\t* Optimal delta-v: \t\t\t\t", delta_v, "m/s")
                    print("\nPropellant Overview of 2nd Stage:")
                    print("\t* Propellant left: \t\t\t\t", m_propellant_left, "kg")
                    print("\t* Propellant used: \t\t\t\t", m_propellant_used, "kg")
                    print("\t* Propellant required by circularization:\t", m_propellant_required, "kg")
                    print("\t* Total propellant used: \t\t\t", m_propellant_total_used_2nd_stage, "kg")
                    print("\t* Possible Payload: \t\t\t\t", (par_roc.M_PROP_2 - m_propellant_total_used_2nd_stage), "kg")

                    # ----- Simulate the rest of the trajectory -----
                    # 1. Coasting
                    # Cutoff second stage engine
                    second_stage_cutoff = True
                    # Define new initial state
                    initial_state_3 = sol_2.y[:, -1]
                    
                    # Define time of simulation 3
                    init_time_3 = sol_2.t[-1]
                    # Calculate time until apogee where we do the delta v burn
                    # print("\n\n")
                    # print("Perigee: ", r_peri_stop)
                    # print("Apogee: ", r_apo_stop)
                    # print("Orbital Period: ", orbit_period_stop)
                    # print("\n\n")
                    time_3 = solvers.get_time_until_apogee(e_stop, initial_state_3[3], initial_state_3[2], orbit_period_stop, a_stop, initial_state_3[1])
                    # print("Time 3: \t\t", time_3)
                    
                    # Call simulation
                    # print("Third Simulation started!")
                    sol_3 = simulate_trajectory(init_time_3, time_3, initial_state_3, False, False, ss_throttle, initial_kick_angle)

                    # 2. Circularization burn
                    initial_state_4 = sol_3.y[:, -1]
                    initial_state_4[2] += delta_v

                    # Calculate time needed to perform the delta v burn
                    burn_time_delta_v = solvers.calculate_burn_time(initial_state_4[4], delta_v)
                    print("\nBurn times:")
                    print("\t* Time for delta-v:\t\t\t\t", burn_time_delta_v, "s\n")

                    initial_state_4[4] -= m_propellant_required

                    # 3. Simulation after circularization burn
                    # Define time of simulation 4
                    init_time_4 = sol_3.t[-1]
                    time_4 = init.DURATION_AFTER_SIMULATION
                    
                    # Call simulation
                    # print("Fourth Simulation started!")
                    sol_4 = simulate_trajectory(init_time_4, time_4, initial_state_4, False, False, ss_throttle, initial_kick_angle)

                    # Collect data and time steps
                    data = np.concatenate((sol_1.y, sol_2.y, sol_3.y, sol_4.y), axis=1)
                    time_steps_simulation = np.concatenate((sol_1.t, sol_2.t, sol_3.t, sol_4.t))

                    # print("Initial time 2: ", init_time_2)
                    # print("End time 2: ", init_time_2 + time_2)
                    # print("Initial time 3: ", init_time_3)
                    # print("End time 3: ", init_time_3 + time_3)
                    # print("Initial time 4: ", init_time_4)
                    # print("End time 4: ", init_time_4 + time_4)
                    # print("\n")
        
                    return time_steps_simulation, data, alt_stop, delta_v, m_propellant_total_used_2nd_stage
   
            else:
                return time_steps_simulation, data, None, 9999999.0, 9999999.0

        elif init.SYM_TYPE == 5:
            return time_steps_simulation, data, double_burn_smallest_prop

        return time_steps_simulation, data
    
    else:
        # Cutoff second stage engine
        second_stage_cutoff = True
        
        #===================================================
        # Simulation after second engine burnout
        #===================================================
        
        # Define new initial state
        initial_state_3 = sol_2.y[:, -1]
        
        # Define time of simulation 3
        init_time_3 = sol_2.t[-1]
        time_3 = init.DURATION_AFTER_SIMULATION
        
        # Call simulation
        # print("Third Simulation started!")
        sol_3 = simulate_trajectory(init_time_3, time_3, initial_state_3, False, False, ss_throttle, initial_kick_angle)
        
        data = np.concatenate((sol_1.y, sol_2.y, sol_3.y), axis=1)
        time_steps_simulation = np.concatenate((sol_1.t, sol_2.t, sol_3.t))
        
        return time_steps_simulation, data
    




def simulate_full_double_burn_trajectory(sol_1, initial_kick_angle):
    """
    Simulates the full double burn trajectory for the rocket after finding the optimal kick angle and optimal time to stop the second engine.

    Input:
        - time_1: time of the first simulation
        - sol_1: solution of the first simulation
        - initial_kick_angle: initial kick angle
        - ss_throttle: throttle of the second stage

    Output:
        - time_steps_simulation: time steps of the full simulation
        - data: data of the full simulation
    """
    global second_stage_cutoff

    # ---- Phase 2: Burn 2nd stage until coasting starts ----
    # Define new initial state
    initial_state_2 = sol_1.y[:, -1]

    # Adjust mass -> perform stage separation
    initial_state_2[4] = initial_state_2[4] - par_roc.M_STRUCTURE_1
    
    # Define time of simulation 2
    init_time_2 = sol_1.t[-1]

    time_2 = coasting_double_burn.double_burn_time_global - init_time_2
    
    # Call simulation for stage 1
    sol_2 = simulate_trajectory(init_time_2, time_2, initial_state_2, False, True, 1.0, initial_kick_angle)


    # ---- Phase 3: Coasting until apogee ----
    par_roc.C_D = 0.

    r_stop = sol_2.y[1, -1]
    v_stop = sol_2.y[2, -1]
    gamma_stop = sol_2.y[3, -1]
    
    # Calculate orbital elements at stop
    a_stop, e_stop, r_apo_stop, r_peri_stop, orbit_period_stop = get_orbital_elements(r_stop, v_stop, gamma_stop)
    
    # Cutoff second stage engine
    second_stage_cutoff = True

    # Define new initial state
    initial_state_3 = sol_2.y[:, -1]
    
    # Define time of simulation 3
    init_time_3 = sol_2.t[-1]

    # Calculate time until apogee where we do the delta v burn
    time_3 = solvers.get_time_until_apogee(e_stop, initial_state_3[3], initial_state_3[2], orbit_period_stop, a_stop, initial_state_3[1])

    # Call simulation
    sol_3 = simulate_trajectory(init_time_3, time_3, initial_state_3, False, False, 1.0, initial_kick_angle)


    # ---- Phase 4: 1st delta v burn ----
    initial_state_4 = sol_3.y[:, -1]
    initial_state_4[2] += coasting_double_burn.double_burn_delta_v_1_global
    m_prop_required_delta_1 = sol_3.y[4, -1] * (1 - np.exp(-coasting_double_burn.double_burn_delta_v_1_global / (c.G0 * par_roc.ISP_2)))

    # Calculate time needed to perform the delta v burn
    burn_time_delta_v_1 = solvers.calculate_burn_time(initial_state_4[4], coasting_double_burn.double_burn_delta_v_1_global)

    initial_state_4[4] -= m_prop_required_delta_1

    # ---- Phase 5: Coasting until transfer apogee ----
    # Define time it takes from the delta v burn until the transfer apogee ---> half of the orbit period
    a_transfer, e_transfer, r_apo_transfer, r_peri_transfer, orbit_period_transfer = get_orbital_elements(initial_state_4[1], initial_state_4[2], initial_state_4[3])
    time_4 = orbit_period_transfer / 2.

    # Define time of simulation 2
    init_time_4 = sol_3.t[-1]
    
    sol_4 = simulate_trajectory(init_time_4, time_4, initial_state_4, False, True, 1.0, initial_kick_angle)


    # ---- Phase 6: 2nd delta v burn ----
    initial_state_5 = sol_4.y[:, -1]
    initial_state_5[2] += coasting_double_burn.double_burn_delta_v_2_global
    m_prop_required_delta_2 = sol_4.y[4, -1] * (1 - np.exp(-coasting_double_burn.double_burn_delta_v_2_global / (c.G0 * par_roc.ISP_2)))

    # Calculate time needed to perform the delta v burn
    burn_time_delta_v_2 = solvers.calculate_burn_time(initial_state_4[4], coasting_double_burn.double_burn_delta_v_2_global)

    initial_state_5[4] -= m_prop_required_delta_2

    # ---- Phase 7: Coasting in final circular orbit ----
    # Define time of simulation 2
    init_time_5 = sol_4.t[-1]
    time_5 = init.DURATION_AFTER_SIMULATION

    # Call simulation for stage 1
    sol_5 = simulate_trajectory(init_time_5, time_5, initial_state_5, False, True, 1.0, initial_kick_angle)


    # ---- Collect data and time steps ----
    data = np.concatenate((sol_1.y, sol_2.y, sol_3.y, sol_4.y, sol_5.y), axis=1)
    time_steps_simulation = np.concatenate((sol_1.t, sol_2.t, sol_3.t, sol_4.t, sol_5.t))

    # ----- Calculate total propellant required -----
    m_propellant_left = sol_2.y[4, -1] - (par_roc.M_STRUCTURE_2 + par_roc.M_PAYLOAD)
    m_propellant_used = par_roc.M_PROP_2 - m_propellant_left
    m_propellant_required = m_prop_required_delta_2 + m_prop_required_delta_1

    # Check if the propellant required is less than the propellant left
    if m_propellant_required < m_propellant_left:
        m_propellant_total_used_2nd_stage = m_propellant_used + m_propellant_required

    print("\nPropellant Overview of 2nd Stage:")
    print("\t* Propellant left: \t\t\t\t", m_propellant_left, "kg")
    print("\t* Propellant used: \t\t\t\t", m_propellant_used, "kg")
    print("\t* Propellant required by circularization:\t", m_propellant_required, "kg")
    print("\t* Total propellant used: \t\t\t", m_propellant_total_used_2nd_stage, "kg")
    print("\t* Possible Payload: \t\t\t\t", (par_roc.M_PROP_2 - m_propellant_total_used_2nd_stage), "kg")

    # print("\n")
    print("\nBurn times:")
    print("\t* Time for delta-v 1:\t\t\t\t", burn_time_delta_v_1, "s")
    print("\t* Time for delta-v 2:\t\t\t\t", burn_time_delta_v_2, "s\n")

    return time_steps_simulation, data

