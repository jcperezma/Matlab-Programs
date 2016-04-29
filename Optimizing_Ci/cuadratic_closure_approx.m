function [a4] = cuadratic_closure_approx(a2)
    ndim=size(a2,1);
    delta=eye(ndim,ndim);
    a4=zeros(ndim,ndim,ndim,ndim);
    for l=1:ndim
        for k=1:ndim
            for j=1:ndim
                for i=1:ndim
                        a4(i,j,k,l)=a2(i,j)*a2(k,l);
                end
            end
        end
    end
end