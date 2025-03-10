""" ===============================================
                Plotting
=============================================== """

import matplotlib.pyplot as plt
import numpy as np

import params.constants as c


def single_run(time_steps, data, angle_attack_values, angle_attack_times):
    """
    Inputs:
        - time_steps: array of time steps (for the data array); [s]
        - data: array of data points. The data array has the following structure:
            * data[0]: downtrack s; [m]
            * data[1]: current radius r from Earth's center; [m]
            * data[2]: velocity norm; [m/s]
            * data[3]: flight path angle; [rad]
            * data[4]: mass of the rocket; [kg]
        - angle_attack_values: array of angle of attack values; [rad]
        - angle_attack_times: array of time steps corresponding to the angle of attack values; [s]

    Plots the following data over time:
        - altitude over downtrack
        - downtrack over time
        - altitude over time
        - velocity norm  over time
        - flight path angle (gamma) over time
        - mass of the rocket over time
        - dynamic pressure over time (based on velocity norm)
        - angle of attack over time
    """

    # -------------- Prepare data --------------
    h = (data[1] - c.r_earth) / 1000.       # altitude h; [km]
    s = data[0] / 1000.                     # downtrack s; [km]

    q = []          # dynamic pressure; [Pa]
    for i in range(len(time_steps)):
        alt = h[i] * 1000
        v = data[2][i]
        rho = c.rho0 * np.exp(-alt / c.H)
        q.append(0.5 * rho * v**2)


    # -------------- Plotting --------------
    fig1, axs1 = plt.subplots(2, 4, figsize=(15, 15))

    # Position plot: r over s
    axs1[0, 0].plot(s, h)
    axs1[0, 0].set_xlabel('downtrack s [km]')
    axs1[0, 0].set_ylabel('altitude h [km]')
    axs1[0, 0].set_title('Trajectory of Rocket')
    axs1[0, 0].grid()

    # Position plot: downtrack over time
    axs1[0, 1].plot(time_steps, s)
    axs1[0, 1].set_xlabel('time [s]')
    axs1[0, 1].set_ylabel('downtrack s [km]')
    axs1[0, 1].set_title('Downtrack over Time')
    axs1[0, 1].grid()

    # Position plot: y over time
    axs1[0, 2].plot(time_steps, h)
    axs1[0, 2].set_xlabel('time [s]')
    axs1[0, 2].set_ylabel('altitude h [km]')
    axs1[0, 2].set_title('Altitude over Time')
    axs1[0, 2].grid()

    # Velocity plot
    axs1[0, 3].plot(time_steps, data[2])
    axs1[0, 3].set_xlabel('time [s]')
    axs1[0, 3].set_ylabel('v [m/s]')
    axs1[0, 3].set_title('Velocity Norm over Time')
    axs1[0, 3].grid()

    # Flight path angle plot
    axs1[1, 0].plot(time_steps, data[3])
    axs1[1, 0].set_xlabel('time [s]')
    axs1[1, 0].set_ylabel('gamma [rad]')
    axs1[1, 0].set_title('Flight Path Angle over Time')
    axs1[1, 0].grid()

    # Mass plot
    axs1[1, 1].plot(time_steps, data[4])
    axs1[1, 1].set_xlabel('time [s]')
    axs1[1, 1].set_ylabel('mass [kg]')
    axs1[1, 1].set_title('Mass of Rocket over Time')
    axs1[1, 1].grid()

    # Dynamic Pressure plot
    axs1[1, 2].plot(time_steps, q)
    axs1[1, 2].set_xlabel('time [s]')
    axs1[1, 2].set_ylabel('q [Pa]')
    axs1[1, 2].set_title('Dynamic Pressure over Time')
    axs1[1, 2].grid()

    # Angle of Attack plot
    axs1[1, 3].plot(angle_attack_times, angle_attack_values)
    axs1[1, 3].set_xlabel('time [s]')
    axs1[1, 3].set_ylabel('angle of attack [deg]')
    axs1[1, 3].set_title('Angle of Attack over Time')
    axs1[1, 3].grid()

    plt.tight_layout()
    plt.show()