# Simulation of a Two Stage Rocket Launch Performing a Gravity-Turn

## Structure of the Project
- **src**: contains the code for the simulation:
  - **components**: contains two files to simulate the environment and the rocket
  - **params**: contains a file for defining constants for the whole project and two files for specifiying parameters for the simulation (one for single run, one for multi run / optimization)
  - **plotting**: contains a file to plot the data for given launch data
  - **simulations**: contains a file that performs the simulation of one launch for given parameters
  - _main\_single\_no\_check.py_: executes a single run without checking final conditions and cretates plots for parameters specified in params.params_single_run.py
  - _main\_single\_with\_check.py_: NOT IMPLEMENTED YET
  - _main\_multi.py_: NOT IMPLEMENTED YET
- **references**: contains references listed below

## Unit Convention
In the code variables are interpreted and used in the respective SI units. The usage of other units only occur when defining input values or plotting the data. This is indicated at the respective part.

## Required Packages
Python packages required to run the simulation:
- numpy
- matplotlib
- scipy

## Notes:
- Instead of determining the altitude where the initial kick takes place, currently a time stamp is defined at which it will be performed. By that the integration can be performed to that time stamp. Otherwise the integration can only be performed for small time steps and the current altitude must be checked with the specified kick altitude.

## To-Do's/Questions:
- final conditions check
- multi_run / optimization
- do we have to include the rotational speed of the Earth that the rocket has during take off? If so, how?
    or does it just reduces the required delta-v for the desired orbit?

## References
Can be found in the [references](/references/) folder:

[1] [Astronautics - The Physics of Space Flight (Ulrich Walter, 2024)](/references/Astronautics%20-%20The%20Physics%20of%20Space%20Flight%20(Ulrich%20Walter,%202024).pdf)