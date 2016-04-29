function [DA_Dt] = Folgar_Tucker_a2_change_dan(W,A, lambda ,D, A4,CI)
    ndim=size(A,1);
    I=eye(ndim);
    gamma_dot_magnitude=sqrt(0.5*(sum(sum(D.^2))));
    %gamma_dot_magnitude=sqrt((sum(sum(D.^2))));
    %disp(gamma_dot_magnitude);
    %DA_Dt=zeros(ndim,ndim);
    [double_dot]=double_dot_product(D,A4);
    DA_Dt= -0.5*(W*A-A*W)+0.5*lambda*(D*A+A*D-2*double_dot)+2*CI*gamma_dot_magnitude*(I-ndim*A);
    
    
end