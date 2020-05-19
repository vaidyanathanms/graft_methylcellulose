# graft_methylcellulose
grafted_methylcellose

source files for grafted methylcellulose systems - FORTRAN/Python

ran_numbers.f90 - for generating random numbers
lammps_inp2.f90 - to create datafile
create_infile.f90 - to generate the MC interaction parameters. Requires temperature as an input. Requires the file cgparams.txt
params.f90 - parameter file for lammps_inp2.f90
main.f90 - Main FORTRAN code for analyzing trajectory files. Requires anainp_mc_var.txt as input

Python files

NOTE: Change "PATH" wherever necessary

genconf3.py - creates initial configuration/runs from restart for different systems. See individual lines for more details. Requires resubmit.py and my_python_functions.py

ana2.py - to analyze LAMMPS trajectories
copy_restart.py/copy_toanalyze.py - to copy restart/analyzed files to system.
NOTE: Other python and FORTRAN files are either deprecated or not used for this proj
