## Implementation of the IOD algorithm

### Overview

This package provides an implementation of the IOD algorithm by Tillman and Spirtes [2011] (http://proceedings.mlr.press/v15/tillman11a/tillman11a.pdf).

The IOD learns equivalence classes of acyclic models  
with latent and selection variables from multiple datasets with overlapping variables. 
It outputs a list of PAGs including the true PAG, if the combined statistics are faithful. 

### Package structure 
The algorithm is implemented in the files IOD.R and IOD_Helper.R. In Simulations_Helper.R, several functions were created to test the 
algorithm and generate simulations. The code has been extensively tested, resulting in numerous files where the results were evaluated, 
leading to subsequent corrections and optimizations of the algorithm. To reproduce the test results from the thesis, the Reproduce_Thesis_results 
folder was created, which describes what was used and how to make the results reproducible.
