% Compute orientation results

a_ij  =  compute_a_ij_fromInitialPositions( 'coords.txt'  );

fprintf(' a_11 = %f  \n',a_ij(1))
fprintf(' a_22 = %f  \n',a_ij(4))
fprintf(' a_33 = %f  \n',a_ij(6))