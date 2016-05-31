
total_time = 5; % total simulation time
dt= 0.1;     % time step

Zeta = 0.1; % parameter
Cb = 0.01;  % parameter
S = 0.5;       % parameter
n_bins = 50;
nFibers = 10000;
l_min = 1e-4; % minimum fiber length [m]
l_max = 15e-3;  %maximum fiber length [m]

Eta = 80; % matrix viscosity [Pa-s]
gamma_dot = 10; % Shear rate [s^-1]
E_f = 73e9;   % Young's modulus of the fiber [Pa]
d_f = 19e-6;  % Diameter of the fiber [m]
p = [Zeta Cb S E_f dt n_bins nFibers l_min l_max Eta gamma_dot total_time d_f];

[Win1,lin1,Ln1,Lw1,cumm_hist1,Ln_hist1,Lw_hist1,t_hist] = fiberAttritionOne(p);