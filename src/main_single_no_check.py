""" ===============================================
     MAIN FILE FOR SINGLE RUN - WITHOUT CHECK
=============================================== """

"""
NOTE:
    - please specify the parameters in the file /parameters/params_single_run.py
    - the simulation will be run without any check of the final conditions
    - plots of the collected data will be created
"""

import params.params_single_run as p
import simulations.single_run as single_run
import plotting.plot_results as plot_results

if __name__ == '__main__':

    # Collect the parameters from the parameters file
    args_single_run = [p.Isp_1, p.F_thrust_1, p.m_structure_1, p.m_prop_1, p.Isp_2, p.F_thrust_2, p.m_structure_2, p.m_prop_2, p.init_angle_of_attack, p.t_initial_kick, p.duration_initial_kick, p.A, p.c_D, p.c_L]

    # Execute the single run simulation
    time_steps, data, angle_attack_times, angle_attack_values = single_run.single_run(args_single_run, False)
    # Execute the plotting of the results
    plot_results.single_run(time_steps, data, angle_attack_values, angle_attack_times)
