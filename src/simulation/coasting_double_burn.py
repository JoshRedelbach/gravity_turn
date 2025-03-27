""" ===============================================
              COASTING DOUBLE BURN
=============================================== """

import components.rocket as rocket
import plotting.plot_results as plot_results
import components.constants as c
import components.solvers as solvers
import numpy as np
import init

def plot(time, data, initial_kick_angle):

    a, e, r_apo, r_peri, _ = rocket.get_orbital_elements(data[1,-1], data[2,-1], data[3,-1])
    
    print("Final orbital elements:")
    print("\t* Semimajor axis:\t", a, "m")
    print("\t* Eccentricity:\t\t", e)
    print("\t* Apoapsis altitude:\t", (r_apo - c.R_EARTH)/1000, "km")
    print("\t* Periapsis altitude:\t", (r_peri - c.R_EARTH)/1000, "km")
    print("\n")

    #===================================================
    # Plot and analyze results
    #===================================================

    # Plot results
    plot_results.single_run(time, data, initial_kick_angle)
    plot_results.plot_trajectory_xy(data)


def execute():

    global DOUBLE_BURN_FULL_SIMULATION, double_burn_smallest_prop_global, double_burn_altitude_global, double_burn_time_global, double_burn_delta_v_1_global, double_burn_delta_v_2_global

    DOUBLE_BURN_FULL_SIMULATION = False
    double_burn_smallest_prop_global = 999999999.
    double_burn_altitude_global = None
    double_burn_time_global = None
    double_burn_delta_v_1_global = None
    double_burn_delta_v_2_global = None

    kick_angle = solvers.find_initial_kick_angle_coast_double_burn()

    print("\nResults:")
    print("\t* Optimal Kick angle: \t\t\t", np.rad2deg(kick_angle), "deg")
    print("\t* Optimal altitude to stop burning:\t", double_burn_altitude_global/1000., "km")
    print("\t* Optimal time to stop burning:\t\t", double_burn_time_global, "s")
    print("\t* Optimal delta-v-total:\t\t", (double_burn_delta_v_1_global + double_burn_delta_v_2_global), "m/s")
    print("\t* Optimal delta-v-1:\t\t\t", double_burn_delta_v_1_global, "m/s")
    print("\t* Optimal delta-v-2:\t\t\t", double_burn_delta_v_2_global, "m/s")

    DOUBLE_BURN_FULL_SIMULATION = True

    time, data = rocket.run(1.0, kick_angle)

    plot(time, data, kick_angle)