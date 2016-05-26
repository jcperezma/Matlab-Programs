#Optimizing Ci

The main program is in OptimizeCIandKappa.m. It uses a Newton Raphson method to fit the RSC model to results of simulations by adjusting Ci and Kappa for each component of the orientation tensor. 

The user has to change manually the following variables: dt, shear_rate, write_freq, max_deformation, delta_time and closureType . The first three variables can be pulled from the fibers.in simulation file, whereas the latter three have to be determined by the user. delta_time is the timestep used in the RSC model, small time steps are required to have a stable solution. There are 10 possible closure approximation types that can be used, it is recommended that the type 9 "ORE" by Cintra is used. The type of closure approximation used has to be reported always with the result. 
