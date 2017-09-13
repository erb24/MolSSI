# MolSSI
Sample codes for MolSSI proposal     

sort.py -- Sorts the modes based on rescaled mode relaxation times.  

barrier_sort.f95 --Sorts LE4PD free-energy surfaces based on timescale of the characteristic motion described by the mode.   

mode_sort.f95 -- Converts the re-ordered modes from reals (given by Python) to integers and writes them to file.   

sort.sh -- Executes all the sorting codes above.   

proj_pca.f95 -- Code to 1) find and diagonalize the covariance matrix from an MD simulation, 2) write the trajectories of these normal modes to file, 3) calculate averages and fluctuations of normal modes, and 4) compare to some properties of the LE4PD.   

PE_generator.nb -- Mathematica notebook to visualize free-energy surfaces and barriers output from the LE4PD.  

contact-map.f95 -- FORTRAN code to generate a contact map for a protein from an MD simulation tracjectory.  

