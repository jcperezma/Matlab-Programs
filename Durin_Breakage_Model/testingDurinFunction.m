clear;
clc;
% Constants
dt = .1; % time step
endTime = 40; % length of simulation desired
L_min = .00001;% min fiber length
L_max = .0005; % max fiber length [m]
dL = L_min; % discrete fiber length intervals [m]
eta = 2000; % viscosity [Pa]
gammaDot = 50; % shear rate [1/s]
E = 72.4*10^9; % Young's Modulus [Pa]
b = 5*10^(-6); % fiber radius [meters]
mExp = 3; % exponent used for P_i, never lower than 3
C_i = .001; % Constant in dPsi_dt function
c{1} = [dt endTime L_min L_max b eta gammaDot E C_i];
% angle range and step size
numAngles = 180;
angles=-pi/2:pi/(numAngles-1):pi/2;
angles = angles';
dphi = pi/numAngles;


% Initial mass length distribution
x = L_min:L_min:L_max;
norm = normpdf(x,.00035,.00015); % first entry is x values, second is mean,
                                 % third is standard deviation
norm = norm/sum(norm);
L_i = norm; 

% Initializing angular orientation
[ psi,phi ] = ComputeFolgarTuckerShearOri( gammaDot, C_i, numAngles,  .00001, 10);
tempPsi1 = psi(1:floor(numAngles/2));
tempPsi2 = psi(floor(numAngles/2)+1:end);
psi = [tempPsi2 tempPsi1]; % to account for how angular distribution used in main program
psi = psi/sum(psi);
c{2} = psi;
p = 3;

% Function Testing
[ L_out ] = DurinFiberBreakage( L_i, p, c );
