""" ===============================================
            SINGLE RUN SIMULATION 
=============================================== """

import params.params_rocket as par_roc
import params.params_simulation as par_sim
import components.rocket as rocket
import plotting.plot_results as plot_results
import params.constants as c
import numpy as np


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
    
    SS_throttle = 2
    initial_kick_angle = - np.deg2rad(17.416805744171135)
    time, data = rocket.run(SS_throttle, initial_kick_angle)
    plot(time, data, par_sim.max_angle_of_attack)
