""" ===============================================
            SINGLE RUN SIMULATION 
=============================================== """

import params.params_rocket as par_roc
import params.params_simulation as par_sim
import components.rocket as rocket
import plotting.plot_results as plot_results
import params.constants as c
import numpy as np


def run(SS_throttle, initial_kick_angle):

    # ---- Debugging ---- 
    # Print desired orbit
    r_desired = c.r_earth + par_sim.alt_desired
    v_desired = np.sqrt(c.mu_earth / r_desired)
    print(r_desired)
    print(v_desired)

    #===================================================
    # Simulation until stage separation
    #===================================================

    # Define initial state
    initial_mass = par_roc.m_structure_1 + par_roc.m_prop_1 + par_roc.m_structure_2 + par_roc.m_prop_2 + par_roc.m_payload
    initial_state_1 = [0., c.r_earth, 0., np.deg2rad(90.), initial_mass]

    # Define time of simulation 1
    time_1 = 500.   #<------TODO

    # Call simulation for stage 1
    sol_1 = rocket.simulate_trajectory(0, time_1, initial_state_1, True, SS_throttle, initial_kick_angle)

    
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
    print("Second Simulation started!")
    sol_2 = rocket.simulate_trajectory(init_time_2, time_2, initial_state_2, False, SS_throttle, initial_kick_angle)
    
    data = np.concatenate((sol_1.y, sol_2.y), axis=1)
    time_steps_simulation = np.concatenate((sol_1.t, sol_2.t))
    
    plot_results.single_run(time_steps_simulation, data, initial_kick_angle)
    
    return time_steps_simulation, data


def plot(time, data, initial_kick_angle):
    a, e, r_apo, r_peri = rocket.get_orbital_elements(data[1,-1], data[2,-1], data[3,-1])
    
    print("\nSemimajor axis: ", a, "m")
    print("Eccentricity: ", e)
    print("Apoapsis altitude: ", (r_apo - c.r_earth)/1000, "km")
    print("Periapsis altitude: ", (r_peri - c.r_earth)/1000, "km")

    #===================================================
    # Plot and analyze results
    #===================================================

    # Plot results
    plot_results.single_run(time, data, initial_kick_angle)



if __name__ == '__main__':
    
    SS_throttle = 1.1
    initial_kick_angle = - np.deg2rad(80)
    time, data = run(SS_throttle, initial_kick_angle)
    plot(time, data, par_sim.max_angle_of_attack)