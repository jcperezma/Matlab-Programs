function [ L_out ] = DurinFiberBreakage( L_i, p, c )
% L_i is the initial mass fraction by length distribution
% p is the parameters to be optimized, mExp is the only such variable in
% this function
% c{1} is all the parameters that are constant
% c{2} is the angular orientation vector, make sure it runs from -pi/2:pi/2
% Constants
dt = c{1}(1); % time step
endTime = c{1}(2); % length of simulation desired
L_min = c{1}(3);% min fiber length
L_max = c{1}(4); % max fiber length [m]
dL = L_min; % discrete fiber length intervals [m]
b = c{1}(5); % fiber radius [meters]
numLengths = floor((L_max-L_min)/dL)+1;
eta = c{1}(6); % viscosity [Pa]
gammaDot = c{1}(7); % shear rate [1/s]
E = c{1}(8); % Young's Modulus [Pa]
mExp = p; % exponent used for P_i, never lower than 3
C_i = c{1}(9); % Constant in dPsi_dt function
% angle range and step size
numAngles = 180;
angles=-pi/2:pi/(numAngles-1):pi/2;
angles = angles';
dphi = pi/numAngles;
% fiber orientation distribution, from input
psi = c{2};
% initial mass distribution matrix, from input
m{1} = L_i;
% calculation of matrix size and initializing matrices
% fiber lengths list
fiberLengths = zeros(numLengths,1);
L = L_min;
for i = 1:numLengths
    fiberLengths(i) = L;
    L = L+dL;
end
M = zeros(numLengths);
tau = zeros(numLengths,1);
P_ij = zeros(numLengths);
P_l = zeros(numLengths);
P_Bu = zeros(numAngles,numLengths);
Bu = zeros(numAngles,numLengths);
Beta = zeros(numLengths,1);
Beta_e = zeros(numLengths,1);


% Calculating values that don't change with time
for i = 1:numLengths
    % variables that change with length
    a = fiberLengths(i)/2; % fiber length
    Beta(i) = a/b; % fiber aspect ratio
    % lambda = (Beta^2-1)/Beta^2;
    Beta_e(i) = Beta(i)*.74;
    
    %%%%%%%%%%%%% Calculating P_Bu for all angle values %%%%%%%%%%%%%%%
    % Creating matrix
    
    M_Ang = zeros(numAngles,1);
    if Beta(i) > 2.88
        for j = 1:numAngles
            M_Ang(j) = sin(angles(j))*cos(angles(j));
            sigma_b = eta*gammaDot*M_Ang(j)*Beta(i)^2/(log(2*Beta(i))-1.75);
            sigma_m = -E*(pi/(4*Beta(i)))^2;
            Bu(j,i) = sigma_b/sigma_m;
            if Bu(j,i) < 0
                Bu(j,i) = 0;
            end
            if Bu(j,i) < 1
                P_Bu(j,i) = (1-exp(-Bu(j,i)))/(1-exp(-1));
            else
                P_Bu(j,i) = 1;
            end
        end
    end
    
    %%%%%%%%% Calculation of P_l(x) for all possible lengths %%%%%%%%%
    if Beta(i) > 2.88
        x = fiberLengths(1:i)-fiberLengths(i)/2;
        for j = 1:i
            P_l(i,j) = 1-exp(-(1-(x(j)/(fiberLengths(i)/2))^2)^mExp);
        end
        totP_l = sum(P_l(i,:));
        P_l(i,:) = P_l(i,:)/totP_l;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%% Calculation of P_tr %%%%%%%%%%%%%%%%%%%%%%%
    P_tr(i) = 0; % initializing for sum
    for j = 1:numAngles
        P_tr(i) = P_tr(i)+ P_Bu(j,i)*psi(j);%*dphi;
    end
    
    %%%%%%%%%%%%%%%%%%%%%% Calculation of tau_j %%%%%%%%%%%%%%%%%%%%%%
    t_r(i) = (2*pi/gammaDot)*(Beta_e(i)+1/Beta_e(i));
    tau(i) = P_tr(i)/t_r(i);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%% Calculation of P_ij %%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:numLengths
    % j refers to fiber length L_j being broken up into fiber length L_i
    for j = 1:numLengths
        P_ij(i,j) = 2*P_l(i,j);
    end
end


% Calculating values that change with time
for k = 1:floor(endTime/dt)
    
    %%%%%%%%%%%%%%%%%%%%%%%%% Calculation M_ij %%%%%%%%%%%%%%%%%%%%%%%%
    for j = 1:numLengths
        M(i,j) = (fiberLengths(i)/fiberLengths(j))*m{k}(j);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% find dm(i) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dm{k} = zeros(numLengths,1);
    for i = 1:numLengths
        dm{k}(i) = 0;
        for j = i:numLengths
            dm{k}(i) = dm{k}(i) + tau(j)*P_ij(j,i)*M(i,j)*dt;
        end
        dm{k}(i) = dm{k}(i)-tau(i)*m{k}(i)*dt;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% update m(i) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 1:numLengths
        m{k+1}(i,1) = m{k}(i)+dm{k}(i);
    end
end
L_out = m{end};
end

