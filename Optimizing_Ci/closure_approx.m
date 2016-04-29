function [a4] = closure_approx(a2)
ndim=size(a2,1);
if ndim==2;
    coeff1=1/24;
    coeff2=1/6;
    A=2;
    B=1;
else
    coeff1=1/35;
    coeff2=1/7;
    A=1.5;
    B=0.5;
end
delta=eye(ndim,ndim);
a4_hat=zeros(ndim,ndim,ndim,ndim);
a4_tilde=a4_hat;
for l=1:ndim
    for k=1:ndim
        for j=1:ndim
            for i=1:ndim
                    a4_hat(i,j,k,l)=(-coeff1)*(delta(i,j)*delta(k,l)...
                                           + delta(i,k)*delta(j,l)...
                                           + delta(i,l)*delta(j,k))...
                                           + (coeff2)*(a2(i,j)*delta(k,l)...
                                                     + a2(i,k)*delta(j,l)...
                                                     + a2(i,l)*delta(j,k)...
                                                     + a2(k,l)*delta(i,j)...
                                                     + a2(j,l)*delta(i,k)...
                                                     + a2(j,k)*delta(i,l));
                    a4_tilde(i,j,k,l)=a2(i,j)*a2(k,l);
            end
        end
    end
end
%disp(a4_hat);

f=0;
 
for j=1:ndim
    for i=1:ndim
         f=f+A*a2(i,j)*a2(i,j);
     end
end

f=f-B;
a4=a4_hat;
%a4=a4_tilde;
%a4=(1-f).*a4_hat+f*a4_tilde;

end