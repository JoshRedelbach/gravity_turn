""" ===============================================
              COASTING SINGLE BURN
=============================================== """

import components.rocket as rocket
import plotting.plot_results as plot_results
import components.constants as c
import components.solvers as solvers
import numpy as np


def plot(time, data, initial_kick_angle):
    a, e, r_apo, r_peri, _ = rocket.get_orbital_elements(data[1,-1], data[2,-1], data[3,-1])
    
    print("Final Orbital Elements")
    print("\t* Semimajor axis:\t\t\t\t", a, "m")
    print("\t* Eccentricity:\t\t\t\t\t", e)
    print("\t* Apoapsis altitude:\t\t\t\t", (r_apo - c.R_EARTH)/1000, "km")
    print("\t* Periapsis altitude:\t\t\t\t", (r_peri - c.R_EARTH)/1000, "km")
    print("")

    #===================================================
    # Plot and analyze results
    #===================================================

    # Plot results
    plot_results.single_run(time, data, initial_kick_angle)
    plot_results.plot_trajectory_xy(data)


def execute():

    global SINGLE_BURN_FULL_SIMULATION, TIME_TO_STOP_BURNING_SINGLE_BURN_FINAL
    
    SINGLE_BURN_FULL_SIMULATION = False
    TIME_TO_STOP_BURNING_SINGLE_BURN_FINAL = None

    kick_angle = solvers.find_initial_kick_angle_coast_single_burn()
    print("\nResults:")
    print("\t* Optimal kick angle: \t\t\t\t", np.rad2deg(kick_angle), "deg")  

    SINGLE_BURN_FULL_SIMULATION = True
    time, data, alt_stopped, delta_v, m_propellant_total_used_2nd_stage = rocket.run(1.0, kick_angle)

    plot(time, data, kick_angle)