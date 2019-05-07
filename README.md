# Value_of_Info_Social_Ecol_Systems
Matlab and R code to replicate analysis from article: "General rules for environmental management to prioritise social-ecological systems research based on a value of information approach."

=======DESCRIPTION=======

This repository contains (1) MATLAB code files to evaluate the expected value of perfect information on "X" (EVPXI), an unknown parameter of interest, in four coupled, social-ecological system; and (2) R code file to recreate Figure 4 from the manuscript. 

The files comprise MATLAB scripts:
1_Code_VOI
2_Code_EVPXI

And one R script:
3_Code_Figure

These files should be run in sequence. 1_Code_VOI generates a series of MATLAB Data files, which are aggregated in 2_Code_EVPXI to calculate EVPXI for each system. 2_Code_EVPXI generates 4 csv files which are combined to reproduce Figure 4 from the manuscript in 3_Code_Figure.

=======REQUIREMENTS=======

MATLAB, preferably version 9.6.0.1072779 (R2019a) or later
R, preferably version 3.52 or later; and R library: grDevices 

Refer any questions to:
Email: k.davis@uq.edu.au
