""" ===============================================
            SINGLE RUN SIMULATION 
    !! NO CHECK OF FINAL CONDITIONS YET !!
=============================================== """

import params.params_rocket as par_roc
import components.rocket as rocket
import plotting.plot_results as plot_results
import params.constants as c
import numpy as np


def run():

    # Define initial state
    initial_mass = par_roc.m_structure_1 + par_roc.m_prop_1 + par_roc.m_structure_2 + par_roc.m_prop_2
    initial_state = [0., c.r_earth, 0., np.deg2rad(90.), initial_mass]

    # Define time of simulation
    time = par_roc.t_burn_1 + par_roc.t_burn_2

    # Execute simulation
    sol = rocket.simulate_trajectory(time, initial_state)
    print("\n\nSimulation finished successfully!\n\n")

    # Plot results
    plot_results.single_run(sol.t, sol.y)



if __name__ == '__main__':
    run()