function [double_dot] = double_dot_product(D,A4)
%Computes the double dot product between a second order tensor and a fourth
ndim=size(D,1);
double_dot=zeros(ndim,ndim);
for j=1:ndim
    for i=1:ndim
        for k=1:ndim
            for l=1:ndim
                double_dot(i,j)=double_dot(i,j)+D(l,k)*A4(i,j,k,l);
            end
        end
    end
end

