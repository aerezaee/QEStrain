# QEStrain
A simple Python program for calculating the mechanical coefficients. 
Just prepare the base input file (it could be scf or relax, that's your call), put the PS files into ps directory, use config.ini to configure run and then run the script by "python StrainCal.py". 
This script makes the structures and inputs base on the prepared input file, put the finished runs into "output" directory (also use these outputs to continue intrupted run) and finally make an final output files and plot the predicted polynomial with calculated errors. 
