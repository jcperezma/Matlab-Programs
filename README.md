# Matlab-Programs

This repo contains:


[Postprocessing Functions](Durin_Breakage_Model)
--------------------
**compute_length_fromFile.m**: Finds the length distributions, average length by number and by weight through time from simulation results

**TrackLoci.m**: Tracks the end of a fiber. Useful to measure rotation periods. 

**compute_a_ij_fromFile.m**: Computes the 3D orientation tendor components through time from simulation results

**compute_a_ij2D_fromFile.m**: Computes the 2D orientation tendor components through time from simulation results, useful to compare results with Folgar and Tucker experimental results

[Optimize Ci and Kappa](Optimizing_Ci)
--------------------
The functions in this folder fit the RSC model to the results obtained with the simulation. the main file is **OptimizeCIandKappa.m** 

[RSC model](RSC model)
--------------------
Functions that computes the development of orientation using the RSC model. The function to call is **ComputeRSC_Evolution**. There is a test function **TestClosures.m** that computes the orientation evolution with different closure approximations and serves as an example of how to call the function. By making Kappa 0 the model is the standard folgar tucker model with Advani's tensorial representation.

[Phelps' Breakage model](Phelps_Breakage_Model)
--------------------
Function that implements the breakage model from Phelps

[Durin's Breakage model](Durin_Breakage_Model)
--------------------
Function that implements the breakage model from Durin

