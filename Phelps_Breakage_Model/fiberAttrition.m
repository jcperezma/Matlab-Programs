function [accum,l_i,Ln,Lw,cumm_hist,Ln_hist,Lw_hist,t_hist] = fiberAttrition(p)

Zeta = p(1); % parameter
Cb = p(2);  % parameter
dt= p(5);     % time step
total_time = p(12); % total simulation time
S = p(3);       % parameter
n_bins = p(6);
nFibers = p(7);
l_min = p(8); % minimum fiber length [m]
l_max = p(9);  %maximum fiber length [m]

Eta_m = p(10); % matrix viscosity [Pa-s]
gamma_dot =p(11); % Shear rate [s^-1]
E_f = p(4);   % Young's modulus of the fiber [Pa]
d_f = p(13);  % Diameter of the fiber [m]


% initialize bins and lengths
l_i = l_min: (l_max-l_min)/(n_bins-1):l_max;

%N_i = ones(1,n_bins)*nFibers/n_bins; % homogeneous distribution
N_i = zeros(1,n_bins);
N_i(n_bins) =nFibers;
Bi =  4*Zeta*Eta_m*gamma_dot*l_i.^4 / (3.14*E_f*d_f^4); % eq 13
P_i = zeros(1,n_bins);
for i = 1: n_bins
   if Bi(i) >= 1
   P_i(i) = Cb * gamma_dot * (1-exp(1-Bi(i)));       % eq 21
   end
end
Lt= sum(N_i.*l_i); % Total fiber length, it has to be preserved during the simulation
Ln = sum( N_i.*l_i )/sum(N_i); % average length by number
Lw = sum( N_i.*l_i.*l_i )/sum(N_i.*l_i); % average length by weight Lw>Ln
% figure(1)
% plot(l_i,N_i,'b')
% hold on
maxN_i = max(N_i);
y = 0: maxN_i/10 :maxN_i;
x_1 = ones(length(y))* Lw;
x_2 = ones(length(y))* Ln;
% plot(x_1, y,'bo' )
% plot(x_2, y,'bx' )
% plot([Lw Lw], [0, max(N_i)],'bo' )
% plot([Ln Ln], [0, max(N_i)],'bx' )
%Compute and normalize R_ik, it is constant through the simulation 

R_ik = zeros(n_bins);

for i = 1:n_bins
   for k=i+1:n_bins
    R_ik(i,k) = normpdf(l_i(i),l_i(k)/2,S*l_i(k)); % eq 22
       
   end
end

alpha_k =zeros(1,n_bins);
for j = 2 :n_bins
   R_sum = sum( R_ik(:,j));
   alpha_k(j) = 2* P_i(j) /R_sum;
end

time = 0;
dNidt =zeros(1,n_bins);
% initialize and compute historic variables.

frame=1;
Lt= sum(N_i.*l_i)
Ln_hist(frame) = sum( N_i.*l_i )/sum(N_i)
Lw_hist(frame) = sum( N_i.*l_i.*l_i )/sum(N_i.*l_i)
cumm_hist(frame,:) =N_i; 
t_hist(frame)=time;
while time<total_time
    for i = 1 : n_bins
        RHS_sum =0;
        for k = i+1 : n_bins
            RHS_sum = RHS_sum + R_ik(i,k) * alpha_k(k) * N_i(k);
        end
        dNidt(i) = -P_i(i)*N_i(i) + RHS_sum;
     end
    N_i = N_i + dt *(dNidt  ); % explicit euler
    time = time + dt;
    if mod1(time,0.1) <0.001
        % compute cummulative length distribution
         cumm  = findCummulative( N_i,l_i );
         cumm_hist(frame,:) = cumm;
        % Coompute length averages
        Lt= sum(N_i.*l_i);
        frame= frame+1;
        Ln_hist(frame) = sum( N_i.*l_i )/sum(N_i);
        Lw_hist(frame) = sum( N_i.*l_i.*l_i )/sum(N_i.*l_i);
        
        %2
        t_hist(frame)=time;
        
        N_i_hist(frame,:) =N_i; 
    end
end

% for i =1 : length(N_i)
%     xmin = N_i(i)*0.6;
%     xmax  = N_i(i)*1.4;
%     x=xmin+rand(1,1)*(xmax-xmin) ;
%     N_i(i) = x;
% end

Wt = sum(N_i.*l_i);
Wi = N_i.*l_i /Wt;

accum(1) = Wi(1);
for i =2 : length(N_i)
    accum(i) = Wi(i) +accum(i-1);
end



% figure(1)
% hold on
% plot(l_i,N_i,'r')
%Lt= sum(N_i.*l_i);
Ln = sum( N_i.*l_i )/sum(N_i);
Lw = sum( N_i.*l_i.*l_i )/sum(N_i.*l_i);
% 
% % plot([Lw Lw], [0, max(N_i)],'ro' )
% % plot([Ln Ln], [0, max(N_i)],'rx' )
% maxN_i = max(N_i);
% y = 0: maxN_i/10 :maxN_i;
% x_1 = ones(length(y))* Lw;
% x_2 = ones(length(y))* Ln;
% plot(x_1, y,'--ro' )
% plot(x_2, y,'--rx' )
