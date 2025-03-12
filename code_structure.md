# Notes for the Structure of the Code

## Two Integrations

1. First engine:
    - until one of the following events occured:
      - current radius larger than radius of desired orbit
      - stage separation done, meaning:
         * first egnine cut-off -> mass propellant 1 is burnt
         * delay after cut-off past
         * ... (maybe now)

Assume instantaneous change of the rocket mass by subtracting the mass of structure 1

1. Second engine:
   - or one of the following events occured:
     - successfully reached desired orbit (check of radius, velocity and flight path angle)
     - no fuel left