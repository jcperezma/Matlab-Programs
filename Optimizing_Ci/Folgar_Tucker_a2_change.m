function [Da2_Dt] = Folgar_Tucker_a2_change(a2, a4, omega, gamma_dot, lambda, C1)
ndim=size(a2,1);

gamma_dot_magnitude=sqrt(0.5*(sum(sum(gamma_dot.^2))));
delta=eye(ndim,ndim);
Da2_Dt=zeros(ndim,ndim);
    accum=zeros(ndim,ndim);
    for j=1:ndim
        for i=1:ndim
            for k=1:ndim
                for l=1:ndim
                    accum(i,j)=accum(i,j)+gamma_dot(k,l)*a4(i,j,k,l);
                end
            end
            Da2_Dt(i,j)=-0.5.*(omega(i,:)*a2(:,j)-a2(i,:)*omega(:,j))...
                        +0.5*lambda*(gamma_dot(i,:)*a2(:,j)+a2(i,:)*gamma_dot(:,j))...
                        -lambda*accum(i,j)...
                        +2*C1*gamma_dot_magnitude*(delta(i,j)-ndim*a2(i,j));
                        
        end
    end
    
end

