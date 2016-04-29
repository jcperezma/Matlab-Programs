function [ error ] = computeOrientationError( phi_mech, psi_mech, numPoints, phi_calc, psi_calc, dphi)

error = 0;

for i = 1: numPoints
   
    %interpolate psi mech
    
    index1 = floor((phi_mech(i)-( -pi/2))/dphi) + 1; 
         
    index2 = index1 +1;
    %linear interpolation
    
    if index2>  numel(psi_calc)
       index1 =  numel(psi_calc) -1;
       index2 = index1 +1;
        
    end
    
    m = (psi_calc(index2) - psi_calc(index1))/dphi;
    b = psi_calc(index1) - m * phi_calc(index1);
    
    psi2compare  = m * phi_mech(i) + b;
    
    error =  error + (psi_mech(i) - psi2compare)^2;
    
end



end

