function [ Da2_Dt] = ComputeRSC_Change( a2, omega, gamma_dot, Ci, kappa, closureType)
%ComputeRSC_Change Computes the material derivative of the second rank
%orientation tensor
%   Dependencies: closure_approx2, find_M_and_L, RSC_a2_change
%   in
%   a2          :  the current orientation tensor, ndim x ndim matrix
%   omega       :  the vorticity tensor, 3x3 matrix
%   gamma_dot   :  the rate of strain tensor, 3x3 matrix
%   Ci          :  the Folgar-Tucker interaction Coefficient, [0.00001,0.1]
%   kappa       :  the Reduced Strain Closure parameter, (0,1], 1 is
%                  Folgar-Tucker original model
%   closureType :  Closure approximation used to approximate the 4th rank
%                  orientation tensor. Based on Tucker and Advani (1999)
%                  1- Isotropic, 2-Linear, 3-Quadratic, 4-Strong flow, 
%                  5-Composite1(Hinch & Leal), 6-Composite2(Hinch & Leal),
%                  7-Hybrid.
%   out
%   Da_Dt       : rate of change of the second order orientation tensor.

[a4] = closure_approx2(a2,closureType);
[M, L] = find_M_and_L(a2,ndim);
[Da2_Dt] = RSC_a2_change( a2, a4, omega, gamma_dot, lambda, kappa, Ci, ndim, M, L ); 


end

