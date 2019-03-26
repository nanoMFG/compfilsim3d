# compfilsim3d
NanoHUB Composite Modeling Project

For [nanohub.org](https://nanohub.org)

## Overview

This project is intended to be a design tool for anyone interested in making conductive composite 3D printing filaments utilizing nanowires. The simulation uses inputs including the additive material properties, size, and amount (vol %), and estimates bulk electrical resistivity based on a representative volume element. The simulation is run through multiple iterations in a Monte Carlo method to produce results from many different possible nanowire arrays within the material. Results from this simulation can be used to guide experimental formulations of new conductive 3D printing filaments.

## Input Details

Nanowire diameter:
* Diameter of the nanowires, in nanometers

Nanowire length:
* Average length of a nanowire, in millimeters.

% Volume of nanowires:
* Percent of the composite's volume that is composed of nanowires.

Alignment (degrees):
* Single value maximum angle of nanowires to the direction of current. For example, 90 degrees would allow the nanowires to be randomly oriented while 0 degrees would force the nanowires to be parallel.

Alignment range (degrees):
* To graph resistivity for multiple alignment degrees, check the box and use the slider to select a range of values for alignment values. The alignment value is the maximum angle of nanowires to the direction of current. For example, 90 degrees would allow the nanowires to be randomly oriented while 0 degrees would force the nanowires to be parallel.

Nanowire material:
* Choose from a list of materials.

## Simulation Details

In this tool, the simulation begins by randomly generating nanowires inside the specified domain until the desired volume percent is reached. The wires are generated with the desired angular alignment as well. The simulation then finds connecting paths from one end of the polymer matrix to the other end and turns that path into a circuit before calling the SPICE circuit solver on it. The solver is quite time consuming depending on the parameters given. Below are some sample parameters and their associated run times per iteration. Note that the order of growth is not linear with respect to volume.

## Tool Output

2 Types:
1. Single value output. If the checkbox is not selected, output will be a single value of resistivity and a geometric visualization of the nanowires.
2. Graph output. If the checkbox is selected, the output will be a graph of resistivity vs. alignment range value. Alignment range to be graphed is specified by the input slider bar.
