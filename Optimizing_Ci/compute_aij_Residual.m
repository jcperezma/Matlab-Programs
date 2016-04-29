function [ residuals ] = compute_aij_Residual( aij_mech, totStrain_mech, aij_calc, totStrain_calc, componentIndex)

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
    b = aij_calc(index1,componentIndex) - m * totStrain_calc(index1,componentIndex);
    
    a_ij2compare  = m * totStrain_mech(i) + b;
    
    residuals(i) =  (aij_mech(i,componentIndex) - a_ij2compare);
    
end


end
