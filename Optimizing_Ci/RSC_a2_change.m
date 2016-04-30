function [ da_dt ] = RSC_a2_change( A, a4, W, D, lambda, kappa, CI, ndim, M, L )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
ndim=size(a4,1);
    I=eye(ndim);
    gamma_dot_magnitude=sqrt(0.5*(sum(sum(D.^2))));% according to wang is 2 instead of 0.5, page 4 of his PhD thesis
    %gamma_dot_magnitude=sqrt((sum(sum(D.^2))));
    %disp(gamma_dot_magnitude);
    da_dt       =      zeros(ndim,ndim);
    RSC_term    =      a4+(1-kappa)* (L - double_dot_4rd_n_4rd( M, a4 ));
    [double_dot]=      double_dot_product(D,RSC_term);
    da_dt       =      -0.5*(W*A-A*W)+ 0.5*lambda*(D*A+A*D-2*double_dot)+ 2*CI*kappa*gamma_dot_magnitude*(I-ndim*A);

end


