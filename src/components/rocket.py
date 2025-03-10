""" ===============================================
    Functions to simulate rocket's behaviour
=============================================== """

import components.environment as env
import params.constants as c

import numpy as np
from scipy.integrate import solve_ivp


class Rocket:
    
    def __init__(self, Isp, mass, F_thrust, c_D, A):
        self.s = 0.                                                 # current downtrack; [m]
        self.r = c.r_earth                                          # current radius from Earth's center; [m]
        self.v = 0.0                                                # current velocity norm; [m/s]
        self.gamma = np.deg2rad(90.0)                               # current flight path angle; [rad]
        self.alpha = 0.0                                            # current angle of attack; [rad]
        self.Isp = Isp                                              # current specific impulse; [s]
        self.mass = mass                                            # current mass; [kg]
        self.F_thrust = F_thrust                                    # norm of current thrust; [N]
        self.c_D = c_D                                              # drag coefficient; [no unit]
        self.c_L = 0.0                                              # lift coefficient; [no unit]
        self.A = A                                                  # cross sectional area; [m^2]


    def set_state(self, state):
        """ 
        Set the new state of the rocket. 
        
        Input:
            - state: new state vector composing of [s, r, v, gamma, mass]
        """
        self.s, self.r, self.v, self.gamma, self.mass = state

    
    def separate_stage(self, mass, Isp, F_thrust):
        """
        Separates the stage of the rocket by setting the mass of the rocket to new value, as well as updating the Isp and the thrust to corresponding values of the new engine.
        
        Input:
            - mass: new mass of the rocket after separation; [kg]
            - Isp: specific impulse of the new engine; [s]
            - F_thrust: thrust of the new engine; [N]
        """
        self.mass = mass
        self.Isp = Isp
        self.F_thrust = F_thrust


    def initial_kick(self, angle_of_attack, c_L):
        """
        Initiates the kick of the rocket by setting the angle of attack to the desired value.
        
        Input:
            - angle_of_attack: desired angle of attack; [rad]
            - c_L: lift coefficient; [no unit]
        """
        self.alpha = angle_of_attack
        self.c_L = c_L


    def end_kick(self):
        """
        Ends the kick of the rocket by resetting the angle of attack and the lift coefficient to zero.
        
        Input:
            - angle_of_attack: desired angle of attack; [rad]
            - c_L: lift coefficient; [no unit]
        """
        self.alpha = 0.0
        self.c_L = 0.0


    def rocket_dynamics(t, state, F_T, Isp, c_D, c_L, A, alpha):
        """
        Simulates the dynamics of the rocket. This function will be integrated by the scipy.solve_ivp function
        
        Input:
            - t: time variable (necessary for solve_ivp function)
            - state: current state vector
            - F_T: norm of force due to thrust; [N]
            - Isp: specific impulse of the engine; [s]
            - c_D: drag coefficient; [no unit]
            - c_L: lift coefficient; [no unit]
            - A: cross sectional area; [m^2]
            - alpha: angle of attack; [rad]

        NOTE: The state vector is defined as follows:
            - s: downtrack; [m]
            - r: radius from Earth's center; [m]
            - v: velocity norm; [m/s]
            - gamma: flight path angle; [rad]
            - m: current mass; [kg]

        Output:
            - derivatives of the state vector
        """
        # Get state components
        s, r, v, gamma, m = state

        # Compute altitude above Earth's surface
        alt = r - c.r_earth                # altitude of the rocket; [m]

        # --- Determine current accelerations and forces ---
        a_grav = env.accel_grav(r)                                  # gravity at the altitude of the rocket in negative radial direction
        F_D = env.drag_force(v, alt, c_D, A)      # drag force norm acting in negative velocity direction
        F_L = env.lift_force(v, alt, c_L, A)      # lift force norm acting in vertical to velocity direction

        # --- Get trigonometric operations of gamma and alpha ---
        c_gamma = np.cos(gamma)
        s_gamma = np.sin(gamma)
        c_alpha = np.cos(alpha)
        s_alpha = np.sin(alpha)

        # --- Compute the derivatives ---

        dsdt = (c.r_earth / r) * v * c_gamma

        # Reference: [1] Equations 6.3.8
        drdt = v * np.sin(gamma)

        # Reference: [1] Equation 6.3.6
        dvdt = (F_T / m) * c_alpha - (F_D / m) - a_grav * s_gamma

        # Catch the case of zero velocity to avoid division by zero
        epsilon = 1e-6
        if v < epsilon:
            dgammadt = 0.
        else:
            # Reference: [1] Equation 6.3.7
            dgammadt = (1./v) * ( (F_T / m) * s_alpha + F_L / m  - (a_grav - (v**2 / r)) * c_gamma )
        
        dmdt = - F_T / (Isp * c.g0)

        return [dsdt, drdt, dvdt, dgammadt, dmdt]


    def simulate_trajectory(self, time_stamp):
        state_init = [self.s, self.r, self.v, self.gamma, self.mass]

        t_span = (0, time_stamp)
        t_eval = np.linspace(0, time_stamp, 1000)
        
        return solve_ivp(Rocket.rocket_dynamics, y0=state_init, t_span=t_span, t_eval=t_eval, args=(self.F_thrust, self.Isp, self.c_D, self.c_L, self.A, self.alpha), max_step=1)
