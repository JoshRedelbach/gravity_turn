""" ===============================================
                Plotting
=============================================== """

import matplotlib.pyplot as plt
import numpy as np

import components.constants as c
import init as par_sim
import components.rocket as r


def single_run(time_steps, data, INITIAL_KICK_ANGLE):
    """
    Inputs:
        - time_steps: array of time steps (for the data array); [s]
        - data: array of data points. The data array has the following structure:
            * data[0]: downtrack s; [m]
            * data[1]: current radius r from Earth's center; [m]
            * data[2]: velocity norm; [m/s]
            * data[3]: flight path angle; [rad]
            * data[4]: mass of the rocket; [kg]

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

    # Reduce data array
    data_reduced = data[:, ::10]
    time_reduced = time_steps[::10]

    # -------------- Prepare data --------------
    h = (data_reduced[1] - c.R_EARTH) / 1000.       # altitude h; [km]
    s = data_reduced[0] / 1000.                     # downtrack s; [km]

    q = []          # dynamic pressure; [Pa]
    for i in range(len(time_reduced)):
        alt = h[i] * 1000
        v = data_reduced[2][i]
        rho = c.RHO0 * np.exp(-alt / c.H)
        q.append(0.5 * rho * v**2)
    
    # Recreate angle of attack values
    # Initialize empty list with length of t
    angle_of_attacks = [0.0] * len(time_reduced)

    for i, t in enumerate(time_reduced):
        if t < r.time_kick_start:
            angle_of_attacks[i] = 0.0
        elif t > (r.time_kick_start + par_sim.DURATION_INITIAL_KICK):
            angle_of_attacks[i] = 0.0
        elif t > (r.time_kick_start + (par_sim.DURATION_INITIAL_KICK / 2.)):
            angle_rate = (t - (r.time_kick_start + r.time_raise)) / (r.time_raise)
            angle_of_attacks[i] = INITIAL_KICK_ANGLE * (1 - angle_rate)
        else:
            angle_rate = (t - r.time_kick_start) / (r.time_raise)
            angle_of_attacks[i] = INITIAL_KICK_ANGLE * angle_rate

    # -------------- Plotting --------------
    fig1, axs1 = plt.subplots(2, 4, figsize=(15, 15))

    # Position plot: r over s
    axs1[0, 0].plot(s, h)
    axs1[0, 0].set_xlabel('downtrack s [km]')
    axs1[0, 0].set_ylabel('altitude h [km]')
    axs1[0, 0].set_title('Trajectory of Rocket')
    axs1[0, 0].grid()

    # Position plot: downtrack over time
    axs1[0, 1].plot(time_reduced, s)
    axs1[0, 1].set_xlabel('time [s]')
    axs1[0, 1].set_ylabel('downtrack s [km]')
    axs1[0, 1].set_title('Downtrack over Time')
    axs1[0, 1].grid()

    # Position plot: y over time
    axs1[0, 2].plot(time_reduced, h)
    axs1[0, 2].set_xlabel('time [s]')
    axs1[0, 2].set_ylabel('altitude h [km]')
    axs1[0, 2].set_title('Altitude over Time')
    axs1[0, 2].grid()

    # Velocity plot
    axs1[0, 3].plot(time_reduced, data_reduced[2])
    axs1[0, 3].set_xlabel('time [s]')
    axs1[0, 3].set_ylabel('v [m/s]')
    axs1[0, 3].set_title('Velocity Norm over Time')
    axs1[0, 3].grid()

    # Flight path angle plot
    axs1[1, 0].plot(time_reduced, np.rad2deg(data_reduced[3]))
    axs1[1, 0].set_xlabel('time [s]')
    axs1[1, 0].set_ylabel('gamma [rad]')
    axs1[1, 0].set_title('Flight Path Angle over Time')
    axs1[1, 0].grid()

    # Mass plot
    axs1[1, 1].plot(time_reduced, data_reduced[4])
    axs1[1, 1].set_xlabel('time [s]')
    axs1[1, 1].set_ylabel('mass [kg]')
    axs1[1, 1].set_title('Mass of Rocket over Time')
    axs1[1, 1].grid()

    # Dynamic Pressure plot
    axs1[1, 2].plot(time_reduced, q)
    axs1[1, 2].set_xlabel('time [s]')
    axs1[1, 2].set_ylabel('q [Pa]')
    axs1[1, 2].set_title('Dynamic Pressure over Time')
    axs1[1, 2].grid()

    # Angle of Attack plot
    axs1[1, 3].plot(time_reduced, np.rad2deg(angle_of_attacks))
    axs1[1, 3].set_xlabel('time [s]')
    axs1[1, 3].set_ylabel('angle of attack [deg]')
    axs1[1, 3].set_title('Angle of Attack over Time')
    axs1[1, 3].grid()

    plt.tight_layout()
    plt.show()
    
    
    
def plot_trajectory_xy(data):
    """
    Plots the rocket trajectory in x-y coordinates with Earth shown as a blue disk.
    
    Inputs:
        - data: array of data points with the following structure:
            * data[0]: downtrack s; [m]
            * data[1]: current radius r from Earth's center; [m]
    """
    # Reduce data_reduced array
    data_reduced = data[:, ::10]

    # -------------- Prepare data --------------
    h = (data_reduced[1] - c.R_EARTH)       # altitude h; [m]
    s = data_reduced[0]                     # downtrack s; [m]
    
    # Convert to cartesian coordinates
    x, y = r.cartesian_coordinates(h, s)
    x = x/1000.
    y = y/1000.

    # Plot setup
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.set_facecolor("black")

    # Plot trajectory
    ax.plot(x, y, color="white", linewidth=1, label="Rocket Trajectory")
    
    # Create Earth representation (circular disk)
    earth_radius_km = c.R_EARTH / 1000.0
    earth = plt.Circle((0, 0), earth_radius_km, color='blue', zorder=1)
    
    # Show Earth
    ax.add_patch(earth)

    # Labels and aesthetics
    ax.set_xlabel("Downtrack Distance (km)", color="white")
    ax.set_ylabel("Altitude (km)", color="white")
    ax.set_title("Rocket Trajectory", color="white")
    ax.tick_params(colors='white')
    ax.grid(color='gray', linestyle='--', linewidth=0.5)
    ax.legend()

    # Set limits to make sure Earth is fully visible
    ax.set_xlim(min(x) - 1200, max(x) + 1200)  # Adjust margins around trajectory
    ax.set_ylim(min(y) - 1200, max(y) + 1200)
    ax.set_aspect('equal')  # Keep aspect ratio realistic

    plt.savefig("rocket_trajectory.jpg", dpi=1000, bbox_inches="tight", pad_inches=0)
    plt.show()
