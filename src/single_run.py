""" ===============================================
            SINGLE RUN SIMULATION 
=============================================== """

import params.params_rocket as par_roc
import components.rocket as rocket
import plotting.plot_results as plot_results
import params.constants as c
import numpy as np


def run():

    #===================================================
    # Simulation until stage separation
    #===================================================

    # Define initial state
    initial_mass = par_roc.m_structure_1 + par_roc.m_prop_1 + par_roc.m_structure_2 + par_roc.m_prop_2 + par_roc.m_payload
    initial_state_1 = [0., c.r_earth, 0., np.deg2rad(90.), initial_mass]

    # Define time of simulation 1
    time_1 = 0.   #<------TODO

    # Call simulation for stage 1
    sol = rocket.simulate_trajectory(time_1, initial_state_1, True)

    
    #===================================================
    # Simulation after stage separation
    #===================================================
    
    # Define new initial state
    initial_state_2 = sol.y[-1]

    # Adjust mass -> perform stage separation
    initial_state_2[4] = initial_state_2[4] - par_roc.m_structure_1
    
    # Define time of simulation 2
    time_2 = 0.   #<------TODO
    
    # Call simulation for stage 1
    sol = rocket.simulate_trajectory(time_2, initial_state_2, False)


    #===================================================
    # Plot and analyze results
    #===================================================
    
    # Plot results
    plot_results.single_run(sol.t, sol.y)



if __name__ == '__main__':
    run()