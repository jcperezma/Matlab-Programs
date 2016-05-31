function [ psi,phi ] = ComputeFolgarTuckerShearOri( gammaDot, C_i, numAngles,  dt, timeEnd)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
phi=-pi/2:pi/(numAngles-1):pi/2;
dphi = pi/numAngles;
psi = ones(numAngles,1)*(1/numAngles);
psi = psi'/ sum(psi'*pi/numAngles );
sum(psi*dphi);

N =numAngles;
psi=(ones(1,N)).*(1/pi);

C4=1./dphi.^2;
C5=1/(2*dphi);
time=0;
while time < timeEnd
        for j=2:N-1;
           dpsi_dphi(j)=   C5*(psi(j+1)-psi(j-1));
           d2psi_dphi2(j)= C4*(psi(j+1)-2.*psi(j)+psi(j-1));
        end
        %calculate the derivatives taking advantage of the periodicity of psi
       
        dpsi_dphi(1)  =(psi(2)-psi(N-1)).*C5;
        dpsi_dphi(N)  = (psi(2)-psi(N-1)).*C5;
        d2psi_dphi2(1)= (psi(2)-2.*psi(1)+psi(N-1)).*C4;
        d2psi_dphi2(N)= (psi(2)-2.*psi(1)+psi(N-1)).*C4;
        
        
%          dpsi_dphi(1)  =0;
%          dpsi_dphi(N)  = 0;
%          d2psi_dphi2(1)= 0;
%          d2psi_dphi2(N)= 0;

%psi_hist(i,:) = psi;
dPsi_dt = dpsi_dphi.*gammaDot.*sin(phi).^2+ psi.*gammaDot.*sin(2*phi)+C_i*gammaDot*d2psi_dphi2;
psi = psi+dPsi_dt*dt;
time = time+dt;


% apply boundary conditions
 %psi(1) =psi(2);
 %psi(N) = psi(N-1);
psi = psi/sum(psi*dphi);
end

end