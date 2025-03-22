""" ===============================================
                GRAVITY TURN SOLVERS
=============================================== """

import numpy as np
from scipy.optimize import bisect, minimize
from components.rocket import run

import components.constants as c
import init as par_sim
import components.rocket as rocket

# =======================================================
#  Direct no coast orbit injection solver (gravity turn)
# =======================================================

def kick_angle_objective(kick_angle, throttle):

    """
    Objective function to find the initial kick angle for the gravity turn.
    
    Input:
        - kick angle: initial kick angle; [rad]
        - throttle: throttle setting for the second stage
    
    Output:
        - last_gamma: flight path angle at the end of the simulation; [rad]
    """
    time, data = run(throttle, kick_angle)
    last_gamma = data[3, -1]
    return last_gamma

def find_initial_kick_angle(throttle):
    """
    Finds the initial kick angle for the gravity turn using a bounded optimization method.
    
    Input:
        - throttle: throttle setting for the second stage
    
    Output:
        - initial kick angle: initial kick angle; [rad]
    """
    bounds = [(-np.deg2rad(90), -np.deg2rad(0.1))]
    result = minimize(lambda x: abs(kick_angle_objective(x[0], throttle)), x0=[-np.deg2rad(10)], bounds=bounds, tol=1e-7)
    return result.x[0]

def second_throttle_objective(throttle):
    print("----------------------------")
    print("Throttle: ", throttle)
    INITIAL_KICK_ANGLE = find_initial_kick_angle(throttle)
    time, data = run(throttle, INITIAL_KICK_ANGLE)
    last_radius = data[1, -1]
    
    # Calculate the maximum altitude reached during the simulation
    max_radius = np.max(data[1, :])
    
    last_gamma = data[3,-1]
    # print("Initial conditions: ", data[:,0])
    # print("Final conditions", data[:,-1])
    print("Kick angle", np.rad2deg(INITIAL_KICK_ANGLE))
    print("Last gamma: ", np.rad2deg(last_gamma))
    print("Throttle: ", throttle)
    a, e, r_apo, r_peri = rocket.get_orbital_elements(data[1, -1], data[2, -1], data[3, -1])
    print("Perigee: ", r_peri - c.R_EARTH, "m")
    print("Apogee: ", r_apo - c.R_EARTH, "m")
    print("Eccentricity: ", e)
    print("Last Altitude: ", (last_radius - c.R_EARTH)/1000, "km")
    print("Max Altitude: ", (max_radius - c.R_EARTH)/1000, "km")
    r_desired = c.R_EARTH + par_sim.ALT_DESIRED
    v_desired = np.sqrt(c.MU_EARTH / r_desired)
    print("delta_Velocity: ", data[2, -1] - v_desired, "m/s")
    
    return max_radius - c.R_EARTH - par_sim.ALT_DESIRED

def find_throttle():
    return bisect(second_throttle_objective, 0.6, 2, xtol=1e-6, maxiter=500)



# =======================================================
#  Hohman Solver
# =======================================================
    
def circularize_delta_v(r, v):
    """
    NOTE: NOT TESTED YET!

    Computes the delta-v required to circularize an orbit at radius r with velocity v.

    Input:
        - r: radius of the orbit; [m]
        - v: velocity of the spacecraft; [m/s]
    
    Output:
        - delta_v: delta-v required to circularize the orbit; [m/s]
    """
    # Compute required velocity at circular orbit with radius r
    v_circular = np.sqrt(c.MU_EARTH / r)

    # Compute delta-v
    delta_v = v_circular - v

    return delta_v


def hohman_transfer(v1, r1, r2):
    """
    NOTE: NOT TESTED YET!

    Computes the delta-v required to perform first burn at r1 to set new apogee to r2 and second burn at r2 to circularize the orbit.

    Input:
        - v1: velocity of the spacecraft at r1; [m/s]
        - r1: radius of the initial orbit; [m]
        - r2: radius of the final orbit; [m]

    Output:
        - delta_v1: delta-v required for the first burn; [m/s]
        - delta_v2: delta-v required for the second burn; [m/s]
        - delta_v_total: total delta-v required for the transfer; [m/s]
    """
    # Compute semi-major axis of the transfer orbit
    a_transfer = (r1 + r2) / 2

    # Compute required velocity of the spacecraft at r2
    v2 = np.sqrt(c.MU_EARTH * (2 / r2 - 1 / a_transfer))

    # Compute delta-v required for the first burn
    delta_v1 = v2 - v1

    # Compute delta-v required for the second burn
    delta_v2 = circularize_delta_v(r2, v2)
    
    # Compute total delta-v required for the transfer
    delta_v_total = delta_v1 + delta_v2

    return delta_v1, delta_v2, delta_v_total