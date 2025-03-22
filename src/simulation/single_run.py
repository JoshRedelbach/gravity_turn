""" ===============================================
            SINGLE RUN SIMULATION 
=============================================== """

import init
import components.rocket as rocket
import plotting.plot_results as plot_results
import components.constants as c
import numpy as np


def plot(time, data, initial_kick_angle):
    a, e, r_apo, r_peri = rocket.get_orbital_elements(data[1,-1], data[2,-1], data[3,-1])
    
    print("\nSemimajor axis: ", a, "m")
    print("Eccentricity: ", e)
    print("Apoapsis altitude: ", (r_apo - c.R_EARTH)/1000, "km")
    print("Periapsis altitude: ", (r_peri - c.R_EARTH)/1000, "km")

    #===================================================
    # Plot and analyze results
    #===================================================

    # Plot results
    plot_results.single_run(time, data, initial_kick_angle)
    plot_results.plot_trajectory_xy(data)

    
def execute():
    time, data = rocket.run(init.SS_THROTTLE, init.INITIAL_KICK_ANGLE)
    plot(time, data, init.INITIAL_KICK_ANGLE)