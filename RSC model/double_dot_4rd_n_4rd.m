function [ result ] = double_dot_4rd_n_4rd( M, a4 )
% Computes the iner product of two 4rd order tensors

ndim=size(a4,1);
result= zeros(ndim,ndim,ndim,ndim);

for i=1:ndim
    for j=1:ndim
        for k=1:ndim
            for l=1:ndim
                for m=1:ndim
                    for n=1:ndim
                        result(i,j,k,l)= result(i,j,k,l) + M(i,j,m,n)*a4(n,m,k,l);
                        
                    end
                end
            end
        end
    end
end
                
            




end

