# Simulation of a Two Stage Rocket Launch Performing a Gravity-Turn With Atmosphere

## Authors
This software project was developed by Joshua Redelbach and Alexandre Pereira.

## Description
This software optimizes the angle of attack with respect to the used propellant for two different scenarios:
1. **Coasting Single Burn**: In this scenario, the rocket performs a gravity turn until a certain altitude, stops burning and starts coasting until the apogee is reached. This apogee corresponds to the desired altitude. When reaching this point, a certain $\Delta v$ is applied to circularize the orbit.
2. **Coasting Double Burn**: In this scenario, the rocket performs a gravity turn until a certain altitude, stops burning and starts coasting until the apogee is reached. This apogee does not have to correspond to the desired altitude. When reaching this point, a certain $\Delta v_1$ is applied to get on a transfer orbit which has its apogee at the desired altitude. There, another $\Delta v_2$ is applied to circularize the orbit.

## Usage
After installing the required python packages (listed below), just specify all the parameters in [_init.py_](/src/init.py) and execute the [_main.py_](/src/main.py)-file.


## Structure of the Project
- [**src**](/src/): Contains the code for the simulation:
  - [**components**](/src/components/): Contains in total five files:
    - [_constants.py_](/src/components/constants.py): Defines all constants required within the whole simulation.
    - [_environment.py_](/src/components/environment.py): Simulates the environment for the whole simulation (gravity, drag, lift).
    - [_rocket_selector.py_](/src/components/rocket_selector.py): Only for handling the different launch vehicle modules. (-> can be ignored from user perspective)
    - [_rocket.py_](/src/components/rocket.py): Implements The main logic of a single launch trajectory simulation, e.g. dynamics of the rockets, integration and interrupt & event handling.
    - [_solvers.py_](/src/components/solvers.py): Implements functions (e.g. cost functions, $\Delta v$ calculations, ...) in order to evaluate and optimize the launch trajectory with respect to different parameters.
  - [**launch_vehicles**](/src/launch_vehicles/): Contains different files which all defines the same parameters and only differ regarding the specific values, which are chosen with respect to different launch vehicles which can be used within the simulation.
  - [**plotting**](/src/plotting/): Contains a file to plot the data for given launch trajectory data.
  - [**simulation**](/src/simulation/): Contains different files, each executing a certain type of simulation.
    - [_single\_run.py_](/src/simulation/single_run.py): Implements a single run for given fixed parameters and plots the results of the launch trajectory.
    - [_direct\_noCoast\_injection.py_](/src/simulation/direct_noCoast_injection.py): Implements the optimization when no coasting phase shall be performed.
  - [_init.py_](/src/init.py): Contains the defintion of all input parameters by the user. For each run of the simulation, the user shall specify the desired parameter values here.
  - [_main.py_](/src/main.py): This file needs to be executed after specifying the parameters and the desired type of simulation specified in _init.py_.
- [**references**](/references/): contains references listed below

## Unit Convention
In the code variables are interpreted and used in the respective SI units. The usage of other units only occur when defining input values or plotting the data. This is indicated at the respective part.

## Required Packages
Python packages required to run the simulation:
- numpy (pip install numpy)
- matplotlib (pip install matplotlib)
- scipy (pip install scipy)

## References
Can be found in the [references](/references/) folder:

[1] [Astronautics - The Physics of Space Flight (Ulrich Walter, 2024)](/references/Astronautics%20-%20The%20Physics%20of%20Space%20Flight%20(Ulrich%20Walter,%202024).pdf)
