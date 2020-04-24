Thermodynamics Project for John Lewis Corker, Somil Joshi, and Chesson 
Sipling. 

This git repo contains code that lets you graph state space with three
variables, two of which are your own choosing. This calculates the 
entropy for the ranges of state variables you give it, and then it also
enables you to view the projection of the state space along with the 
percentage difference between the three different model types. 

The code calculates the entropy for a Van der Waal gas and an Ideal
Gas, however it has two different methods for calculating the entropy
for an Ideal Gas. One way is done macroscopically with thermodynamics, 
while the other is done microscopically with statistical mechanics.

#RUNNING THE CODE#
To generate the plots, you must download the thermo_main_program.py as 
well as graph_class.py and have numpy and matplotlib modules installed
in python. Once this has been done, run the thermo_main_program from 
your command line with the following format:
\n
python3 thermo_main_program.py [State Variable 1] [start of range]
[end of range] [State Variable 2] [start of range] [end of range]
\n
With this done, the plots should automatically be displayed.
