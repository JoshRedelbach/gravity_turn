""" ===============================================
    Functions to simulate the environment
=============================================== """

import numpy as np
import components.constants as c

# -------------- Functions --------------

def drag_force(v, alt, C_D, A, RHO0=c.RHO0, H=c.H):
    """
    Function to get the norm of the atmospheric drag at certain height and speed.
    
    Reference: [1] Equation 6.3.10

    Input:
        - v: norm current velocity vector; [m/s]
        - alt: current altitude above Earth's surface; [m]
        - C_D: drag coefficient; [no unit]
        - A: cross sectional area; [m^2]
        - RHO0: sea-level air density; [kg/m^3]
        - H: scale height of atmosphere; [m]

    Outputs:
        - drag: norm of drag force; [N]
    """
    
    rho = RHO0 * np.exp(-alt / H)
    drag = 0.5 * C_D * rho * v**2 * A
    return drag


def lift_force(v, alt, C_L, A, RHO0=c.RHO0, H=c.H):
    """
    Function to get the norm of the lift at certain height and speed.
    
    Reference: [1] Equation 6.3.9

    Input:
        - v: norm current velocity vector; [m/s]
        - alt: current altitude above Earth's surface; [m]
        - C_L: lift coefficient; [no unit]
        - A: cross sectional area; [m^2]
        - RHO0: sea-level air density; [kg/m^3]
        - H: scale height of atmosphere; [m]

    Outputs:
        - lift: norm of lift force; [N]
    """
    
    rho = RHO0 * np.exp(-alt / H)
    drag = 0.5 * C_L * rho * v**2 * A
    return drag


def accel_grav(r, mu=c.MU_EARTH):
    """
    Function to get the acceleration due to gravity at certain height.
    
    Reference: [1] Equation 6.3.11

    Input:
        - r: current radius of rocket to Earth's center; [m]
        - mu: grav. constant Earth; [m^3/s^2]
    Outputs:
        - a_grav: norm of gravity acceleration; [m/s^2]
    """

    a_grav = mu / ((r)**2)
    return a_grav
