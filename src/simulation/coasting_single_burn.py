""" ===============================================
              COASTING SINGLE BURN
=============================================== """

import components.rocket as rocket
import plotting.plot_results as plot_results
import components.constants as c
import init
import components.solvers as solvers
import numpy as np


def plot(time, data, initial_kick_angle):
    a, e, r_apo, r_peri = rocket.get_orbital_elements(data[1,-1], data[2,-1], data[3,-1])
    
    print("Semimajor axis:\t\t ", a, "m")
    print("Eccentricity:\t\t ", e)
    print("Apoapsis altitude: \t", (r_apo - c.R_EARTH)/1000, "km")
    print("Periapsis altitude: \t", (r_peri - c.R_EARTH)/1000, "km")
    print("\n\n")

    #===================================================
    # Plot and analyze results
    #===================================================

    # Plot results
    # plot_results.single_run(time, data, initial_kick_angle)
    # plot_results.plot_trajectory_xy(data)


def execute():
    # time, data, alt_stopped, delta_v, m_propellant_total_used_2nd_stage = rocket.run(init.SS_THROTTLE, init.INITIAL_KICK_ANGLE)
    # print("Altitude stopped: \t\t", str(alt_stopped/1000), "km")
    # print("Delta V: \t\t\t", str(delta_v/1000), "km/s")
    # print("\n")
    # plot(time, data, init.INITIAL_KICK_ANGLE)

    kick_angle = solvers.find_initial_kick_angle_coast_single_burn()
    print("Kick angle: ", np.rad2deg(kick_angle))    