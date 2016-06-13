function [ a2, a2_in_time, strain ] = ComputeRSC_Evolution( a2, omega, gamma_dot, totalStrain, Ci, kappa, closureType,delta_time )
%ComputeRSC_Evolution Computes the evolution of the second rank
%orientation tensor
%   Dependencies: closure_approx2, find_M_and_L, RSC_a2_change
%   in
%   a2          :  the initial orientation tensor, ndim x ndim matrix
%   omega       :  the vorticity tensor, 3x3 matrix
%   gamma_dot   :  the rate of strain tensor, 3x3 matrix
%   totalStrain :  total deformation the system will undergo
%   Ci          :  the Folgar-Tucker interaction Coefficient, [0.00001,0.1]
%   kappa       :  the Reduced Strain Closure parameter, (0,1], 1 is
%                  Folgar-Tucker original model
%   closureType :  Closure approximation used to approximate the 4th rank
%                  orientation tensor. Based on Tucker and Advani (1999)
%                  1-Isotropic, 
%                  2-Linear, 
%                  3-Quadratic, 
%                  4-Strong flow, 
%                  5-Composite1(Hinch & Leal), 
%                  6-Composite2(Hinch & Leal),
%                  7-Hybrid.
%                  8-ORF.
%                  9-ORE
%                  10-ORW
%   out
%   a2           :  The final orientation tensor
%   a2_in_time   :  the history of a2 
%   strain       :  a vector that holds the values of the strains


gamma_dot_magnitude=sqrt(0.5*(sum(sum(gamma_dot.^2))));
total_time=totalStrain/gamma_dot_magnitude;
n_time_steps=round(total_time/delta_time);
ndim=length(a2(:,1));
lambda =1;
% initialize variables
t=zeros(n_time_steps,1);
a2_in_time=zeros(n_time_steps,ndim^2);

k=1;
    for i=1:ndim
        for j=1:ndim
            a2_in_time(1,k)=a2(i,j);
            k=k+1;
        end
    end
    
for it=2:n_time_steps
    
    t(it,1)=(it-1)*delta_time;
    
    % Compute the 4th order tensor a4, L and M
    [a4] = closure_approx2(a2,closureType);
    [M, L] = find_M_and_L(a2,ndim);
    
    % Compute the rate of change of a2
    [Da2_Dt] = RSC_a2_change( a2, a4, omega, gamma_dot, lambda, kappa, Ci, ndim, M, L );
    
    
    % Update a2 with explicit Euler
    a2(1:ndim,1:ndim)=a2(1:ndim,1:ndim)+delta_time.*Da2_Dt(1:ndim,1:ndim);
    
    % Normalize the diagonal of the tensor
    a2=a2/trace(a2);
        
    k=1;
    for i=1:ndim
        for j=1:ndim
            a2_in_time(it,k)=a2(i,j);
            k=k+1;
        end
    end
    
end
strain=gamma_dot_magnitude.*t(:,1);


end

