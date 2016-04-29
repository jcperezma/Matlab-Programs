
clc
clear all
close all

% Script that finds RSC constants from a mechanistic model Simulation

% STEP 1%
% Load mechanistic model data and transform it into tensors.
folder_name = uigetdir; % select the folder that contains the output folder
numFramesFileName=[folder_name '\output\nbr_frames.txt'];
positionsFileName=[folder_name '\output\positions.out'];
% the data below was pulled from Fibers.in
dt = 1e-5;
shear_rate = 10;
write_freq = 2000;
fprintf('Loading Mechanistic model results \n' );
[a_ij_mech, totalDeformation_mech ] =  compute_a_ij_fromFile(numFramesFileName, positionsFileName, dt, shear_rate, write_freq  );
fprintf('Mechanistic model results loaded \n' );


% STEP 2%
% Initialize values
%%%%%%%%%%%%%%%%%%%%%%%   INPUT   %%%%%%%%%%%%%%%%%%%%%%%
shear_rate =1;  %Shear Flow
epsilon_dot=+1; %Elongational rate
delta_time=0.1;
totalStrain=max(totalDeformation_mech*1.1);
ndim=3;
minDeformation = 50;
maxDeformation = 800;
if ndim == 3
a_ij_calc_indices = [1 2 3 5 6 9]; 
titles = ['a_{11}'; 'a_{12}'; 'a_{13}'; 'a_{22}'; 'a_{23}'; 'a_{33}' ];
end
if ndim == 2
a_ij_calc_indices = [1 2 4]; 
titles = ['a_{11}'; 'a_{12}'; 'a_{22}'];
end
re=20;
Ci=0.1;
flow_type=1;
lambda=1;%(1-re^2)/(1+re^2)
kappa=1;
closureType = 7;
% 3D Random initial orientation
a2=zeros(ndim,ndim);

if ndim == 2 
   a2(1,1)= a_ij_mech(1,1); %a11
   a2(1,2)= a_ij_mech(1,2); %a12
   a2(2,1)= a_ij_mech(1,2); %a12
   a2(2,2)= a_ij_mech(1,3); %a11
end

if ndim ==3
   a2(1,1)= a_ij_mech(1,1); %a11
   a2(1,2)= a_ij_mech(1,2); %a12
   a2(2,1) = a2(1,2);
   a2(1,3)= a_ij_mech(1,3); %a12
   a2(3,1) = a2(1,3);
   a2(2,2)= a_ij_mech(1,4); %a22   
   a2(2,3)= a_ij_mech(1,5); %a23
   a2(3,2) = a2(2,3);
   a2(3,3)= a_ij_mech(1,6); %a33
end

minIndex =1;
for i = 1 : length(totalDeformation_mech)
   if totalDeformation_mech(i) > minDeformation
       minIndex = i;
       break
   end    
end

maxIndex =1;
for i = 1 : length(totalDeformation_mech)
   if totalDeformation_mech(i) > maxDeformation
       maxIndex = i;
       break
   end    
end


[ omega, gamma_dot ] = InitializeFlowVars(shear_rate, epsilon_dot, flow_type,ndim );
%a2=zeros(ndim,ndim);

% Optimization Options
tol = 0.01; % maximum error willing to accept
maxIter =10;
plotIntermidiateResults =0;% 1 if want to plot intermediate results, 0 if not
% Start Loop

% Optimization for a non linear square problem
% you can use Gauss-Newton of Levenberg-Marquardt

min_s = 1000000000000000;
bestBeta= [0 0 ];

numComponents = 3* (ndim-1);

% STEP 3%
% Find Best Ci for each tensorial Component

figure(1)

for component = 1:numComponents
    
    % Number of Parameters of your function
    m = 1; % Ci and Kappa
    % Define inital guess for parameters and deltaBeta
    %            Ci    Kappa
    dBeta    = [0.00001   0.01]; % make this way smaller than the guesses you make
    Beta_old = [0.01    1 ]; %initial guess
    
    
    
    % Set up limit values for the parameters
    %            Ci    Kappa
    Beta_min = [  0.0001      0.01];
    Beta_max = [0.4      1   ];
    
    
    
    useLM = 0; % 0 is Gauss Newton
    % 1 is Levenberg marquardt
    lambda = 0.01;% Initial lambda for LM
    tolLM = 0.0000001;% criterion to accept LM step
    
    % optimization will stop when tolerance or maximum number of iterations are
    % reached
    
    % Data used in the optimization
    % the function is such that y = f(x,Beta)
    % we know y for certain x, we want to find Beta
    % y is y_observed
    
    % y_observed comes from the mechanistic model data, that is a_ij_mech
    y_observed = a_ij_mech(:,component)';
    x = totalDeformation_mech;
    %y_observed = log(yStar);
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%   END OF INPUT   %%%%%%%%%%%%%%%%%%%%%%%
    
    
    % Plot points we want to fit
    %figure(component)
    subplot(ndim,ndim,a_ij_calc_indices(component))
    plot(x,y_observed,'r-X')
    hold on
    
    % Function call!!!!!!!!!
    % you should be able to calculate the value of the function for the x given
    % above and Beta.
    % it is here where you would put your specific function.
    
    %-----------------------------------------------------------                                                    Ci      Kappa
    [ a2_final, a_ij_calc, totalDeformation_calc ] = ComputeRSC_Evolution( a2, omega, gamma_dot, totalStrain, Beta_old(1), Beta_old(2), closureType,delta_time );
    
    x_calculated = totalDeformation_calc;
    
    y_calculated = interpolate_aij_values( a_ij_mech, totalDeformation_mech, a_ij_calc, totalDeformation_calc, a_ij_calc_indices(component));
    %y_calculated =log( Beta_old(1) ) - (Ustar./(R .* (x-Too ))) - (Beta_old(2)./(x.*dT.*f)) ;
    %-----------------------------------------------------------
    
    % Plot first guess
    plot(x,y_calculated,'g' )
    
    % we want to minimize s, the squared error
    %-----------------------------------------------------------
    %s_old = compute_aij_Error( a_ij_mech, totalDeformation_mech, a_ij_calc, totalDeformation_calc, component);
    s_old = sum((y_observed(minIndex:maxIndex)-y_calculated(minIndex:maxIndex)).^2); % squared error
    %-----------------------------------------------------------
    
    clear s_history
    s_count = 1;
    s_history(s_count)=s_old;
    s_count = s_count +1;
    iter = 0;
    min_s = s_old;
    bestBeta =  Beta_old;
    
    while s_old>tol && iter< maxIter
        iter = iter+1;
        Beta_old
        
        % While the error is higher than the tolerance correct parameter values
        
        % STEP 1
        % Find Jacobian
        for j = 1 : m
            Beta_up = Beta_old;
            Beta_up(j) = Beta_up(j) + dBeta(j);
            
            %-----------------------------------------------------------            % Function call!!!!!!!!!
            [ a2_final, a_ij_calc_up, totalDeformation_calc ] = ComputeRSC_Evolution( a2, omega, gamma_dot, totalStrain, Beta_up(1), Beta_up(2), closureType,delta_time );
            y_calculatedUp = interpolate_aij_values( a_ij_mech, totalDeformation_mech, a_ij_calc_up, totalDeformation_calc, a_ij_calc_indices(component));
            %y_calculatedUp =log( Beta_up(1) ) - (Ustar./(R .* (x-Too ))) - (Beta_up(2)./(x.*dT.*f)) ;
            %-----------------------------------------------------------
                        
            % I use forward differences to calculate derivative
            %r_up = compute_aij_Residual( a_ij_mech, totalDeformation_mech, a_ij_calc, totalDeformation_calc, component)/dBeta(j);
            r_up = (y_calculatedUp(minIndex:maxIndex)-y_calculated(minIndex:maxIndex))/dBeta(j);
            
            J(j,:) =  r_up;
        end
        
        % The residuals for current value of parameters.
        
        r = y_observed(minIndex:maxIndex)-y_calculated(minIndex:maxIndex);
        
        
        % STEP 2
        % calculate correction, and find new Beta
        if useLM == 1
            %levenberg-Marquardt correction
            correction = pinv(J*J'+lambda*diag(diag(J*J')))*(J*r');
        else
            % Gauss-Newton Correction
            correction = pinv(J*J')*(J*r');
        end
        
        %update parameter values
        Beta(1:m) = Beta_old(1:m) +  correction';
        Beta(m+1:length(Beta_old)) = Beta_old(m+1:end);
        Beta = min(max(Beta_min,Beta),Beta_max);
        
        % STEP 3
        % calculate function  and error, with new parameters
        %-----------------------------------------------------------        % Function call!!!!!!!!!
        [ a2_final, a_ij_calc, totalDeformation_calc ] = ComputeRSC_Evolution( a2, omega, gamma_dot, totalStrain, Beta(1), Beta(2), closureType,delta_time );
        y_calculated = interpolate_aij_values( a_ij_mech, totalDeformation_mech, a_ij_calc, totalDeformation_calc, a_ij_calc_indices(component));
        %y_calculated =log( Beta(1) ) - (Ustar./(R .* (x-Too ))) - (Beta(2)./(x.*dT.*f)) ;
        %-----------------------------------------------------------
        
        %y_calculated = Beta(1) * cos(Beta(2)*x ) + Beta(2) * sin(Beta(1)*x );
        %y_calculated = Beta(1)*x+Beta(2);
        
        %plot(x,y_calculated )
        
        % Calculate error
        s = sum((y_observed(minIndex:maxIndex)-y_calculated(minIndex:maxIndex)).^2);
        
        % STEP 4
        % if Gauss Newton is used, accept new Beta
        
        if useLM  ~= 1
            Beta_old = Beta;
            s_old = s;
        end
        
        
        % if LM is used, check if step is accepted
        rho = (s_old - s) / ( 2*correction' * (lambda * correction + J*r') ); % Nielsen
        
        if useLM == 1
            if rho > tolLM
                % accept try
                Beta_old = Beta;
                s_old = s;
                fprintf('accept try')
                % decrease lambda
                lambda = max(lambda/9,1.e-7);
            else
                fprintf('do not accept try\n')
                % Increase lambda
                lambda = min(lambda*11,1.e7);
            end
        end
        
        if plotIntermidiateResults ==1
            plot(x,y_calculated );
        end
        
         if s_old < min_s
            min_s = s_old;
            bestBeta =  Beta_old;
        end
        s_history(s_count)=s_old;
        s_count = s_count +1;
        
    end
    fprintf('------------------------------------------------')
    Beta_old = bestBeta;
    [ a2_final, a_ij_calc, totalDeformation_calc ] = ComputeRSC_Evolution( a2, omega, gamma_dot, totalStrain, Beta_old(1), Beta_old(2), closureType,delta_time );
    y_calculated = interpolate_aij_values( a_ij_mech, totalDeformation_mech, a_ij_calc, totalDeformation_calc, a_ij_calc_indices(component));
    %y_calculated =log( Beta_old(1) ) - (Ustar./(R .* (x-Too ))) - (Beta_old(2)./(x.*dT.*f)) ;
    s_old = sum((y_observed(minIndex:maxIndex)-y_calculated(minIndex:maxIndex)).^2); % squared error
    
    fprintf('Component: %d \n', component )
    fprintf('Lowest error is: %f \n', s_old )
    fprintf('Number of iteration steps: %d \n', iter)
    fprintf('Ci: %f \n',Beta_old(1))
    fprintf('Kappa: %f \n ',Beta_old(2))
    fprintf('Error History\n')
    s_history
    subplot(ndim,ndim,a_ij_calc_indices(component))
    plot(x,y_calculated );
    
    str = sprintf(' %s, Best Ci = %f',titles(component,:),Beta_old(1));
    title(str);
    ylim([0 1 ]);
    hold off
    
    % Store the Ci
    BestCis(component) = Beta_old(1);
end
fprintf('Cis found \n')
fprintf('Finding Kappas \n')
% STEP 5
% Find the best Kappa for the Cis found previously.
tol = tol * 0.1; % maximum error willing to accept
saveas(gcf,[folder_name '\output\Optimized Ci.fig'])

figure(2)


% I will swap Ci and Kappa
for component = 1:numComponents
    
    % Number of Parameters of your function
    m = 1; % Ci and Kappa
    % Define inital guess for parameters and deltaBeta
    %            Kappa    Ci 
    dBeta    = [ 0.0001   0.00001]; % make this way smaller than the guesses you make
    Beta_old = [ 0.999      BestCis(component) ]; %initial guess
    
    
    
    % Set up limit values for the parameters
    %            Kappa    Ci
    Beta_min = [  0.01    0.001];
    Beta_max = [  1      0.4];
    
    
    
    useLM = 0; % 0 is Gauss Newton
    % 1 is Levenberg marquardt
    lambda = 0.01;% Initial lambda for LM
    tolLM = 0.0000001;% criterion to accept LM step
    
    % optimization will stop when tolerance or maximum number of iterations are
    % reached
    
    % Data used in the optimization
    % the function is such that y = f(x,Beta)
    % we know y for certain x, we want to find Beta
    % y is y_observed
    
    % y_observed comes from the mechanistic model data, that is a_ij_mech
    y_observed = a_ij_mech(:,component)';
    x = totalDeformation_mech;
    %y_observed = log(yStar);
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%   END OF INPUT   %%%%%%%%%%%%%%%%%%%%%%%
    
    
    % Plot points we want to fit
    %figure(component)
    subplot(ndim,ndim,a_ij_calc_indices(component))
    plot(x,y_observed,'r-X')
    hold on
    
    % Function call!!!!!!!!!
    % you should be able to calculate the value of the function for the x given
    % above and Beta.
    % it is here where you would put your specific function.
    
    %-----------------------------------------------------------                                                    Ci      Kappa
    [ a2_final, a_ij_calc, totalDeformation_calc ] = ComputeRSC_Evolution( a2, omega, gamma_dot, totalStrain, Beta_old(2), Beta_old(1), closureType,delta_time );
    
    x_calculated = totalDeformation_calc;
    
    y_calculated = interpolate_aij_values( a_ij_mech, totalDeformation_mech, a_ij_calc, totalDeformation_calc, a_ij_calc_indices(component));
    %y_calculated =log( Beta_old(1) ) - (Ustar./(R .* (x-Too ))) - (Beta_old(2)./(x.*dT.*f)) ;
    %-----------------------------------------------------------
    
    % Plot first guess
    plot(x,y_calculated,'g' )
    
    % we want to minimize s, the squared error
    %-----------------------------------------------------------
    %s_old = compute_aij_Error( a_ij_mech, totalDeformation_mech, a_ij_calc, totalDeformation_calc, component);
    s_old = sum((y_observed(1:minIndex)-y_calculated(1:minIndex)).^2); % squared error
    %-----------------------------------------------------------
    
    clear s_history 
    clear J
    s_count = 1;
    s_history(s_count)=s_old;
    s_count = s_count +1;
    iter = 0;
    min_s = s_old;
    bestBeta =  Beta_old;
    
    while s_old>tol && iter< maxIter
        iter = iter+1;
        Beta_old;
        
        % While the error is higher than the tolerance correct parameter values
        
        % STEP 1
        % Find Jacobian
        for j = 1 : m
            Beta_up = Beta_old;
            Beta_up(j) = Beta_up(j) + dBeta(j);
            
            %-----------------------------------------------------------            % Function call!!!!!!!!!
            [ a2_final, a_ij_calc_up, totalDeformation_calc ] = ComputeRSC_Evolution( a2, omega, gamma_dot, totalStrain, Beta_up(2), Beta_up(1), closureType,delta_time );
            y_calculatedUp = interpolate_aij_values( a_ij_mech, totalDeformation_mech, a_ij_calc_up, totalDeformation_calc, a_ij_calc_indices(component));
            %y_calculatedUp =log( Beta_up(1) ) - (Ustar./(R .* (x-Too ))) - (Beta_up(2)./(x.*dT.*f)) ;
            %-----------------------------------------------------------
                        
            % I use forward differences to calculate derivative
            %r_up = compute_aij_Residual( a_ij_mech, totalDeformation_mech, a_ij_calc, totalDeformation_calc, component)/dBeta(j);
            r_up = (y_calculatedUp(1:minIndex)-y_calculated(1:minIndex))/dBeta(j);
            
            J(j,:) =  r_up;
        end
        
        % The residuals for current value of parameters.
        
        r = y_observed(1:minIndex)-y_calculated(1:minIndex);
        
        
        % STEP 2
        % calculate correction, and find new Beta
        if useLM == 1
            %levenberg-Marquardt correction
            correction = pinv(J*J'+lambda*diag(diag(J*J')))*(J*r');
        else
            % Gauss-Newton Correction
            correction = pinv(J*J')*(J*r');
        end
        
        %update parameter values
        Beta(1:m) = Beta_old(1:m) +  correction';
        Beta(m+1:end) = Beta_old(m+1:end);
        Beta = min(max(Beta_min,Beta),Beta_max);
        
        % STEP 3
        % calculate function  and error, with new parameters
        %-----------------------------------------------------------        % Function call!!!!!!!!!
        [ a2_final, a_ij_calc, totalDeformation_calc ] = ComputeRSC_Evolution( a2, omega, gamma_dot, totalStrain, Beta(2), Beta(1), closureType,delta_time );
        y_calculated = interpolate_aij_values( a_ij_mech, totalDeformation_mech, a_ij_calc, totalDeformation_calc, a_ij_calc_indices(component));
        %y_calculated =log( Beta(1) ) - (Ustar./(R .* (x-Too ))) - (Beta(2)./(x.*dT.*f)) ;
        %-----------------------------------------------------------
        
        %y_calculated = Beta(1) * cos(Beta(2)*x ) + Beta(2) * sin(Beta(1)*x );
        %y_calculated = Beta(1)*x+Beta(2);
        
        %plot(x,y_calculated )
        
        % Calculate error
        s = sum((y_observed(1:minIndex)-y_calculated(1:minIndex)).^2);
        
        % STEP 4
        % if Gauss Newton is used, accept new Beta
        
        if useLM  ~= 1
            Beta_old = Beta;
            s_old = s;
        end
        
        
        % if LM is used, check if step is accepted
        rho = (s_old - s) / ( 2*correction' * (lambda * correction + J*r') ); % Nielsen
        
        if useLM == 1
            if rho > tolLM
                % accept try
                Beta_old = Beta;
                s_old = s;
                fprintf('accept try')
                % decrease lambda
                lambda = max(lambda/9,1.e-7);
            else
                fprintf('do not accept try\n')
                % Increase lambda
                lambda = min(lambda*11,1.e7);
            end
        end
        
        if plotIntermidiateResults ==1
            plot(x,y_calculated );
        end
        
         if s_old < min_s
            min_s = s_old;
            bestBeta =  Beta_old;
        end
        s_history(s_count)=s_old;
        s_count = s_count +1;
        
    end
    fprintf('------------------------------------------------')
    Beta_old = bestBeta;
    [ a2_final, a_ij_calc, totalDeformation_calc ] = ComputeRSC_Evolution( a2, omega, gamma_dot, totalStrain, Beta_old(2), Beta_old(1), closureType,delta_time );
    y_calculated = interpolate_aij_values( a_ij_mech, totalDeformation_mech, a_ij_calc, totalDeformation_calc, a_ij_calc_indices(component));
    %y_calculated =log( Beta_old(1) ) - (Ustar./(R .* (x-Too ))) - (Beta_old(2)./(x.*dT.*f)) ;
    s_old = sum((y_observed(1:minIndex)-y_calculated(1:minIndex)).^2); % squared error
    
    fprintf('Component: %d \n', component )
    fprintf('Lowest error is: %f \n', s_old )
    fprintf('Number of iteration steps: %d \n', iter)
    fprintf('Ci: %f \n',Beta_old(2))
    fprintf('Kappa: %f \n ',Beta_old(1))
    fprintf('Error History\n')
    s_history
    subplot(ndim,ndim,a_ij_calc_indices(component))
    plot(x,y_calculated );
    
    str = sprintf(' %s, Best Kappa = %f',titles(component,:),Beta_old(1));
    title(str);
    ylim([0 1 ]);
    hold off
    
    % Store the Ci
    BestKappas(component) = Beta_old(1);
end
saveas(gcf,[folder_name '\output\Optimized Kappa.fig'])


% store results in a file


fileID = fopen([folder_name '\output\Best Parameters.txt'],'w');
fprintf(fileID,'%12s %12s %12s\n','Component', 'Best Ci', 'Best Kappa');
for i=1:numComponents
fprintf(fileID,'%12s %12f %12f\n',titles(i,:), BestCis(i) ,BestKappas(i));
end
fclose(fileID);


