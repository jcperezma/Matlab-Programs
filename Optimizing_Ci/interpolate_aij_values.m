function [ aij_interpolated ] = interpolate_aij_values( aij_mech, totStrain_mech, aij_calc, totStrain_calc, componentIndex)

numPoints = length(totStrain_mech);
dt = totStrain_calc(2) - totStrain_calc(1);

for i = 1: numPoints
   
    %interpolate psi mech
    
    index1 = floor((totStrain_mech(i))/dt) + 1; 
         
    index2 = index1 +1;
    %linear interpolation
    
    if index2>  numel(totStrain_calc)
       index1 =  numel(totStrain_calc) -1;
       index2 = index1 +1;
        
    end
    
    m = (aij_calc(index2,componentIndex) - aij_calc(index1,componentIndex))/dt;
    b = aij_calc(index1,componentIndex) - m * totStrain_calc(index1);
    
    aij_interpolated(i)  = m * totStrain_mech(i) + b;
       
end


end