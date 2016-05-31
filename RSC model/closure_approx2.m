function [a4] = closure_approx2(a,closure_approx_kind)
% computes the closure approximations found on Tucker and Advani's paper

ndim=size(a,1);
delta=eye(ndim);

f=1-27*det(a);
cum=0;
for i=1:ndim
    for j=1:ndim
        cum=cum+a(i,j)*a(j,i);
    end
end
alpha=exp((2*(1-3*cum))/(1-cum));


for i=1:8
    coeff(i)=0;
end

switch closure_approx_kind
    case(1)
        coeff(1)=1/15;
        coeff(2)=1/15;
    case(2)
        coeff(1)=-1/35;
        coeff(2)=-1/35;
        coeff(3)=1/7;
        coeff(4)=1/7;
    case(3)
        coeff(5)=1;
    case(4)
        coeff(5)=1;
        coeff(6)=1;
        for p=1:ndim
            for o=1:ndim
                coeff(7)=coeff(7)+(a(o,p)*a(p,o));
               
            end
        end
        coeff(7)=-2/coeff(7);
    case(5)
        coeff(3)=2/5;
        coeff(5)=-1/5;
        coeff(6)=3/5;
        coeff(8)=-2/5;
    case(6)
        
       coeff(1)=26*alpha/315;
       coeff(2)=26*alpha/315;
       coeff(3)=16*alpha/63;
       coeff(4)=-(4*alpha)/21;
       coeff(5)=1;
       coeff(6)=1;
       for p=1:ndim
            for o=1:ndim
                coeff(7)=coeff(7)+(a(o,p)*a(p,o));
               
            end
       end
       coeff(7)=-2/coeff(7);
    case(7)
       coeff(1)=-(1-f)/35;
       coeff(2)=-(1-f)/35;
       coeff(3)=(1-f)/7;
       coeff(4)=(1-f)/7;
       coeff(5)=f;
    case(8)
       a4  = ORF( a );
       return
    case(9)
       a4  = ORE( a );
       return
    case(10)
       a4  = ORW( a );
       return
    case(11)
       a4  = ORW3( a );
       return
end

a4=zeros(3,3,3,3);

for i=1:ndim
    for j=1:ndim
        for k=1:ndim
            for l=1:ndim
                acum1=zeros(ndim,ndim,ndim,ndim);
                acum2=zeros(ndim,ndim,ndim,ndim);
                for m=1:ndim
                    for n=1:ndim
                        acum1(i,j,k,l)=acum1(i,j,k,l)+a(i,m)*a(m,j)*a(k,n)*a(n,l);
                    end
                    acum2(i,j,k,l)=acum2(i,j,k,l)+delta(i,j)*a(k,m)*a(m,l)+a(i,m)*a(m,j)*delta(k,l);
                end
                
                a4(i,j,k,l)=+ coeff(1)* delta(i,j)*delta(k,l)...
                            + coeff(2)*(delta(i,k)*delta(j,l)+delta(i,l)*delta(j,k))...
                            + coeff(3)*(delta(i,j)*a(k,l)+a(i,j)*delta(k,l))...
                            + coeff(4)*(a(i,k)*delta(j,l)+a(j,l)*delta(i,k)+a(i,l)*delta(j,k)+a(j,k)*delta(i,l))...
                            + coeff(5)* a(i,j)*a(k,l)...
                            + coeff(6)*(a(i,k)*a(j,l)+a(i,l)*a(j,k))...
                            + coeff(7)* acum1(i,j,k,l)...
                            + coeff(8)*acum2(i,j,k,l);
            end
         end
    end
end
