close all


re=20;
Ci=0.001;
flow_type=1;
lambda=1;%(1-re^2)/(1+re^2)
kappa=1;

shear_rate =1;  %Shear Flow
epsilon_dot=+1; %Elongational rate
delta_time=0.1;
totalStrain=80;
ndim =3;
a2 = eye(3)*1/3;   
         

figure
hold on
% Hybrid
closureType = 7;
[ omega, gamma_dot ] = InitializeFlowVars(shear_rate, epsilon_dot, flow_type,ndim );
[ a2, a2_in_time, strain] = ComputeRSC_Evolution( a2, omega, gamma_dot, totalStrain, Ci, kappa, closureType,delta_time );
plot(strain, a2_in_time(:,1),'r')

% Composite1(Hinch & Leal)
closureType = 5;
a2 = eye(3)*1/3;
[ omega, gamma_dot ] = InitializeFlowVars(shear_rate, epsilon_dot, flow_type,ndim );
[ a2, a2_in_time, strain] = ComputeRSC_Evolution( a2, omega, gamma_dot, totalStrain, Ci, kappa, closureType,delta_time );
plot(strain, a2_in_time(:,1),'k')

% Composite2(Hinch & Leal)
closureType = 6;
a2 = eye(3)*1/3;
[ omega, gamma_dot ] = InitializeFlowVars(shear_rate, epsilon_dot, flow_type,ndim );
[ a2, a2_in_time, strain] = ComputeRSC_Evolution( a2, omega, gamma_dot, totalStrain, Ci, kappa, closureType,delta_time );
plot(strain, a2_in_time(:,1))

% ORF
closureType = 8;
a2 = eye(3)*1/3;
[ omega, gamma_dot ] = InitializeFlowVars(shear_rate, epsilon_dot, flow_type,ndim );
[ a2, a2_in_time, strain] = ComputeRSC_Evolution( a2, omega, gamma_dot, totalStrain, Ci, kappa, closureType,delta_time );
plot(strain, a2_in_time(:,1),'g')

% ORE
closureType = 9;
a2 = eye(3)*1/3;
[ omega, gamma_dot ] = InitializeFlowVars(shear_rate, epsilon_dot, flow_type,ndim );
[ a2, a2_in_time, strain] = ComputeRSC_Evolution( a2, omega, gamma_dot, totalStrain, Ci, kappa, closureType,delta_time );
plot(strain, a2_in_time(:,1),'--')

% ORW
closureType = 10;
a2 = eye(3)*1/3;
[ omega, gamma_dot ] = InitializeFlowVars(shear_rate, epsilon_dot, flow_type,ndim );
[ a2, a2_in_time, strain] = ComputeRSC_Evolution( a2, omega, gamma_dot, totalStrain, Ci, kappa, closureType,delta_time );
plot(strain, a2_in_time(:,1),'--r')

% ORW3
closureType = 11;
a2 = eye(3)*1/3;
[ omega, gamma_dot ] = InitializeFlowVars(shear_rate, epsilon_dot, flow_type,ndim );
[ a2, a2_in_time, strain] = ComputeRSC_Evolution( a2, omega, gamma_dot, totalStrain, Ci, kappa, closureType,delta_time );
plot(strain, a2_in_time(:,1),'-.')


legend('HYB','HL1','HL2','ORF','ORE','ORW','ORW3')

ylim([0 1])


