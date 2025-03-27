""" ===============================================
                GRAVITY TURN SOLVERS
=============================================== """

import numpy as np
from scipy.optimize import bisect, minimize
from scipy.optimize import differential_evolution
from scipy.optimize import shgo, brute
from scipy.optimize import dual_annealing
from components.rocket import run
import time

import components.constants as c
import init as par_sim
import components.rocket as rocket

import init

from components.rocket_selector import select_rocket
# Selected Launch Vehicle
par_roc = select_rocket(init.LV)

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
    bounds = [(par_sim.ALPHA_LOWEST, par_sim.ALPHA_HIGHEST)]
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
    a, e, r_apo, r_peri, _ = rocket.get_orbital_elements(data[1, -1], data[2, -1], data[3, -1])
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
    delta_v = abs(v_circular - v)

    return delta_v


def hohman_transfer(v1, r1, r2):
    """
    NOTE: NOT TESTED YET!

    Computes the delta-v required to perform first burn at r1 to set new apogee to r2 and second burn at r2 to circularize the orbit.

    Input:
        - v1: velocity of the spacecraft at r1; [m/s]
        - r1: radius of first burn; [m]
        - r2: radius of the final orbit; [m]

    Output:
        - delta_v1: delta-v required for the first burn; [m/s]
        - delta_v2: delta-v required for the second burn; [m/s]
        - delta_v_total: total delta-v required for the transfer; [m/s]
    """
    # Compute semi-major axis of the transfer orbit
    a_transfer = (r1 + r2) / 2.

    # Compute required velocity of the spacecraft at r2
    v2 = np.sqrt(c.MU_EARTH * ((2. / r1) - (1. / a_transfer)))

    # Compute delta-v required for the first burn
    delta_v1 = abs(v2 - v1)

    # Compute velcoity of the spacecraft at r2
    v3 = np.sqrt(c.MU_EARTH * ((2. / r2) - (1. / a_transfer)))

    # Compute delta-v required for the second burn
    delta_v2 = circularize_delta_v(r2, v3)
    
    # Compute total delta-v required for the transfer
    delta_v_total = delta_v1 + delta_v2

    return delta_v_total, delta_v1, delta_v2


def get_time_until_apogee(e, gamma, v, T, a, r):
    """
    NOTE: NOT TESTED YET!

    Computes the time until the spacecraft reaches the apogee of the orbit.

    Input:
        - e: eccentricity of the orbit
        - gamma: flight path angle of the spacecraft; [rad]
        - v: velocity of the spacecraft; [m/s]
        - T: period of the orbit; [s]
        - a: semi-major axis of the orbit; [m]
        - r: radius of the orbit; [m]

    Output:
        - time_until_apogee: time until the spacecraft reaches the apogee of the orbit; [s]
    """
    # print("Radius: ", r)
    # print("Gamma: ", np.rad2deg(gamma))
    # print("Velocity: ", v)

    theta = np.arccos((a * (1 - e**2) - r) / (e * r))                # true anomaly of stop point
    ecc_anomaly = 2 * np.arctan2(np.sqrt((1 - e) / (1 + e)) * (1 - np.cos(theta)), np.sin(theta))    # eccentric anomaly of stop point
    mean_anomaly = ecc_anomaly - e * np.sin(ecc_anomaly)             # mean anomaly of stop point
    time_until_apogee = T / (2 * np.pi) * mean_anomaly               # time until the spacecraft reaches the apogee of the orbit
    time_until_apogee = (T / 2.) - time_until_apogee

    # print("Time until apogee: ", time_until_apogee)

    return time_until_apogee
    

# =======================================================
#  Coasting and Single Burn Solver (gravity turn)
# =======================================================

def coasting_single_burn_objective(kick_angle):
    """
    Objective function to find the initial kick angle for the gravity turn which minimizes the used propellant in the second stage.
    
    Input:
        - kick angle: initial kick angle; [rad]
    
    Output:
        - m_propellant_total_used_2nd_stage: total mass of propellant used in the second stage; [kg]
    """
    time, data, alt_stopped, delta_v, m_propellant_total_used_2nd_stage = run(1.0, kick_angle)

        # ---- Debugging ----
    print("Kick angle:\t\t", np.rad2deg(kick_angle))
    print("Propellant used:\t", m_propellant_total_used_2nd_stage, "kg")
    print("Delta V:\t\t", delta_v, "m/s")
    print("\n")

    # ---- Debugging ----
    # print("\nAltitude stopped: ", str(alt_stopped), "m")
    # print("Delta V: ", str(delta_v/1000), "km/s")
    # print("Kick angle: ", np.rad2deg(kick_angle))
    # print("Propellant used: ", m_propellant_total_used_2nd_stage, "kg")
    # print("\n")

    # a, e, r_apo, r_peri, _ = rocket.get_orbital_elements(data[1,-1], data[2,-1], data[3,-1])
    # print("Semimajor axis:\t\t ", a, "m")
    # print("Eccentricity:\t\t ", e)
    # print("Apoapsis altitude: \t", (r_apo - c.R_EARTH)/1000, "km")
    # print("Periapsis altitude: \t", (r_peri - c.R_EARTH)/1000, "km")
    # print("\n\n")
    return m_propellant_total_used_2nd_stage


def find_initial_kick_angle_coast_single_burn():
    """
    Finds the initial kick angle for the gravity turn using the differential evolution optimizer.
    
    Output:
        - initial kick angle: initial kick angle; [rad]
    """
    bounds = [(par_sim.ALPHA_LOWEST, par_sim.ALPHA_HIGHEST)]

    initial_guess = init.ALPHA_INITIAL_GUESS

    print("\nFinding initial kick angle for coasting single burn...\n")

    # Time measurement
    start_time = time.time()

    # -- DIFFERENTIAL EVOLUTION --
    # result = differential_evolution(
    #     lambda x: abs(coasting_single_burn_objective(x[0])),
    #     bounds=bounds,
    #     tol=1e-7,
    #     strategy='best1bin',
    #     maxiter=1000,
    #     popsize=15,
    #     mutation=(0.5, 1),
    #     recombination=0.7,
    #     x0=[initial_guess]
    # )

    # -- DUAL ANNEALING --
    # result = dual_annealing(
    #     lambda x: abs(coasting_single_burn_objective(x[0])),
    #     bounds=bounds,
    #     maxiter=1000,
    #     initial_temp=5230,  # Default temperature; can be adjusted
    #     visit=2.62,  # Controls exploration (higher = more exploration)
    #     accept=-5.0,  # Controls acceptance probability
    #     x0=[initial_guess]
    # )

    # -- BRUTE FORCE --
    result = brute(
        lambda x: abs(coasting_single_burn_objective(x[0])),
        ranges=bounds,  # Instead of bounds, we use ranges
        Ns=1000,  # Number of grid points per parameter
        finish=None,  # Avoids additional local optimization
        full_output=True
    )
    # print("Result: ", np.rad2deg(result[0]))
    # print("Propellant: ", np.rad2deg(result[1]))


    # -- SHGO --
    # # Use SHGO optimizer
    # result = shgo(
    #     lambda x: abs(coasting_single_burn_objective(x[0])),
    #     bounds=bounds,
    #     sampling_method='sobol'
    # )
    # # Print all local minima found
    # print("All local minima found:")
    # print(result)
    # for i, local_min in enumerate(result.xl):
    #     print(f"Local minimum {i + 1}: Kick angle = {local_min[0]}, Objective value = {result.funl[i]}")


    # Time measurement
    end_time = time.time()
    print(f"Optimization finished after {np.round(end_time - start_time, 2)} seconds.")

    # Return the global minimum
    return result[0]
    return result.x[0]




# =======================================================
#  Coasting and Double Burn Solver (gravity turn)
# =======================================================

def coasting_double_burn_objective(kick_angle):
    """
    Objective function to find the initial kick angle for the gravity turn which minimizes the used propellant in the second stage.
    
    Input:
        - kick angle: initial kick angle; [rad]
    
    Output:
        - double_burn_smallest_prop: total mass of propellant used in the second stage; [kg]
    """
    time_steps_simulation, data, double_burn_smallest_prop = run(1.0, kick_angle)

    # ---- Debugging ----
    print("Kick angle:\t\t", np.rad2deg(kick_angle))
    print("Propellant used:\t", double_burn_smallest_prop, "kg")
    print("\n")

    # a, e, r_apo, r_peri, _ = rocket.get_orbital_elements(data[1,-1], data[2,-1], data[3,-1])
    # print("Semimajor axis:\t\t ", a, "m")
    # print("Eccentricity:\t\t ", e)
    # print("Apoapsis altitude: \t", (r_apo - c.R_EARTH)/1000, "km")
    # print("Periapsis altitude: \t", (r_peri - c.R_EARTH)/1000, "km")
    # print("\n\n")
    return double_burn_smallest_prop



def find_initial_kick_angle_coast_double_burn():
    """
    Finds the initial kick angle for the gravity turn using the differential evolution optimizer.
    
    Output:
        - initial kick angle: initial kick angle; [rad]
    """
    bounds = [(par_sim.ALPHA_LOWEST, par_sim.ALPHA_HIGHEST)]

    initial_guess = init.ALPHA_INITIAL_GUESS

    print("\nFinding initial kick angle for coasting double burn...\n")

    # Time measurement
    start_time = time.time()

    # result = differential_evolution(
    #     lambda x: abs(coasting_double_burn_objective(x[0])),
    #     bounds=bounds,
    #     tol=1e-7,
    #     strategy='best1bin',
    #     maxiter=1000,
    #     popsize=15,
    #     mutation=(0.5, 1),
    #     recombination=0.7,
    #     x0=[initial_guess]
    # )

    # print("Result angle: ", np.rad2deg(result.x[0]))
    # print("Result cost: ", np.rad2deg(result.fun))

    # -- BRUTE FORCE --
    result = brute(
        lambda x: abs(coasting_double_burn_objective(x[0])),
        ranges=bounds,  # Instead of bounds, we use ranges
        Ns=1000,  # Number of grid points per parameter
        finish=None,  # Avoids additional local optimization
        full_output=True
    )

    # Time measurement
    end_time = time.time()
    print(f"Optimization finished after {np.round(end_time - start_time, 2)} seconds.")

    return result[0]

    return result.x[0]




def calculate_burn_time(m0, delta_v):
    """
    Calculates the time the second stage would need in reality to burn the propellant required for a specific delta v.

    Input:
        - delta_v: delta-v required for the burn; [m/s]
        - m0: initial mass of the second stage; [kg]
    
    Output:
        - t_burn: time the second stage would need to burn the propellant; [s]
    """
    t_burn = m0 * ( 1 - np.exp(-delta_v / par_roc.ISP_2 / c.G0) ) / (par_roc.F_THRUST_2 / c.G0 / par_roc.ISP_2)

    return t_burn