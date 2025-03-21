# Simulation of a Two Stage Rocket Launch Performing a Gravity-Turn With Atmosphere

## Authors
This project was developed by Joshua Redelbach and Alexandre Pereira.

## Description
**TO-DO !!**


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
- **references**: contains references listed below

## Unit Convention
In the code variables are interpreted and used in the respective SI units. The usage of other units only occur when defining input values or plotting the data. This is indicated at the respective part.

## Required Packages
Python packages required to run the simulation:
- numpy
- matplotlib
- scipy

## References
Can be found in the [references](/references/) folder:

[1] [Astronautics - The Physics of Space Flight (Ulrich Walter, 2024)](/references/Astronautics%20-%20The%20Physics%20of%20Space%20Flight%20(Ulrich%20Walter,%202024).pdf)
