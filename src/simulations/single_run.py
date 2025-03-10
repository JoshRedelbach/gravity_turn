""" ===============================================
            SINGLE RUN SIMULATION 
    !! NO CHECK OF FINAL CONDITIONS !!
=============================================== """

import components.rocket as r
import params.constants as c

import numpy as np


def single_run(args, check_flag=False):

    """
    Inputs:
        - args: list of arguments
            * args[0]: Isp of the first stage; [s]
            * args[1]: thrust of the first stage; [N]
            * args[2]: mass of the first stage structure; [kg]
            * args[3]: mass of the first stage propellant; [kg]
            * args[4]: Isp of the second stage; [s]
            * args[5]: thrust of the second stage; [N]
            * args[6]: mass of the second stage structure; [kg]
            * args[7]: mass of the second stage propellant; [kg]
            * args[8]: angle of attack for initial kick; [deg]
            * args[9]: time to start the initial kick; [s]
            * args[10]: duration of the initial kick; [s]
            * args[11]: cross sectional area; [m^2]
            * args[12]: drag coefficient; [no unit]
            * args[13]: lift coefficient; [no unit]
            
            When the check of the final conditions is implemented:
            * args[14]: time step for the final conditions check; [s]
            * args[15]: altitude of desired orbit; [m]
            ...
        - check_flag: flag to check the final conditions; [bool]

    Outputs:
        - data: array of the data of the simulation
        - time_steps: array of the time steps of the simulation
        - angle_attack_values: array of the angle of attack values
        - angle_attack_times: array of the time steps of the angle of attack values

    Structure of this simulation:
    
        1. Initialize the rocket with the values of the first stage
        2. Simulate the vertical ascent of the rocket
        3. Simulate the initial kick
        4. Simulate the gravity turn until first stage is burnt out
        5. Separate the first stage
        6. Initialize the rocket with the values of the second stage
        7. Simulate the rest of the gravity turn until second stage is burnt out
           (or desired orbit is reached <--- not implemented yet)
        8. Returns the data and time steps of the simulation

    """
    # Extract arguments
    Isp_1, F_thrust_1, m_structure_1, m_prop_1, Isp_2, F_thrust_2, m_structure_2, m_prop_2, init_angle_of_attack, t_initial_kick, duration_initial_kick, A, c_D, c_L = args

    # Calculate the burn times of the stages
    t_burn_1 = m_prop_1 / (F_thrust_1 / (Isp_1 * c.g0))
    t_burn_2 = m_prop_2 / (F_thrust_2 / (Isp_2 * c.g0))

    # Define array to store the data
    data = None
    time_steps = None

    # 1. Initialize the rocket with the values of the first stage
    m_initial = m_structure_1 + m_prop_1 + m_prop_2 + m_structure_2
    rocky = r.Rocket(Isp_1, m_initial, F_thrust_1, c_D, A)

    # 2. Simulate the vertical ascent of the rocket
    if t_initial_kick > t_burn_1:
        print("Gravity turn has to start before the first stage burnout!")
        exit()
    sol = rocky.simulate_trajectory(t_initial_kick)
    data = sol.y
    # Set the new current state of the rocket to the end of the vertical ascent
    state = [sol.y[0, -1], sol.y[1, -1], sol.y[2, -1], sol.y[3, -1], sol.y[4, -1]]
    time_steps = sol.t
    last_time = sol.t[-1]
    rocky.set_state(state)
    print("\nVertical ascent simulated.")

    # 3. Simulate the initial kick
    rocky.initial_kick(np.deg2rad(init_angle_of_attack), c_L)
    sol = rocky.simulate_trajectory(duration_initial_kick)
    state = [sol.y[0, -1], sol.y[1, -1], sol.y[2, -1], sol.y[3, -1], sol.y[4, -1]]
    rocky.set_state(state)
    data = np.concatenate((data, sol.y), axis=1)
    time_steps = np.concatenate((time_steps, sol.t + last_time))
    last_time = time_steps[-1]
    rocky.end_kick()
    print("\nInitial kick simulated.")

    # 4. Simulate the gravity turn until first stage is burnt out
    sol = rocky.simulate_trajectory(t_burn_1 - (t_initial_kick + duration_initial_kick))
    data = np.concatenate((data, sol.y), axis=1)
    time_steps = np.concatenate((time_steps, sol.t + last_time))
    last_time = time_steps[-1]
    # Set the new current state of the rocket to the end of the first stage burnout
    state = [sol.y[0, -1], sol.y[1, -1], sol.y[2, -1], sol.y[3, -1], sol.y[4, -1]]
    rocky.set_state(state)
    print("\nFirst stage burnt out.")

    # 5. Separate the first stage
    rocky.separate_stage(m_structure_2 + m_prop_2, Isp_2, F_thrust_2)
    print("\nFirst stage separated.")

    # 7. Simulate the rest of the gravity turn until second stage is burnt out or desired orbit is reached
    sol = rocky.simulate_trajectory(t_burn_2)
    data = np.concatenate((data, sol.y), axis=1)
    time_steps = np.concatenate((time_steps, sol.t + last_time))
    state = [sol.y[0, -1], sol.y[1, -1], sol.y[2, -1], sol.y[3, -1], sol.y[4, -1]]
    rocky.set_state(state)
    print("\nSecond stage burnt out.")

    # Create angle of attack values and corresponding time steps
    angle_attack_values = [0, 0, init_angle_of_attack, init_angle_of_attack, 0, 0.]
    angle_attack_times = [0, t_initial_kick, t_initial_kick, t_initial_kick + duration_initial_kick, t_initial_kick + duration_initial_kick, t_burn_1 + time_steps[-1]]
    
    print("\nSimulation completed\n.")

    # Return the data of current run    
    return time_steps, data, angle_attack_times, angle_attack_values