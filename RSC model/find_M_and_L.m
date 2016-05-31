function [ M , L ] = find_M_and_L( a2, n_dim )
%find_M_and_L Returns the fourth-ordered tensors M and L
%   Detailed explanation goes here
   %p - numbers of eigenvalues
%~~~~~~~~~~~~Delete this later~~~~~~~~~~~~~~~~%
% n_dim=3;
% a2=zeros(n_dim,n_dim);
% a2(1,1)=0.5;
% a2(2,2)=0.5;
% a3(3,3)=0.5;

%~~~~~~~~~testing~~~~~~~~~~~~~~~~~%
  
L = zeros(n_dim,n_dim,n_dim,n_dim); %storage term for fourth order tensor L
M = zeros(n_dim,n_dim,n_dim,n_dim); %storage term for fourth order tensor M

[D,V] = eig(a2); %returns matrix D of eigenvalues and matrix V whose 
                %columns are the corresponding right eigenvectors
stor_sigma = diag(V); %returns a column vector of diagonal elements of eigenvalue matrix

   
for i = 1:n_dim
    for j = 1:n_dim
        for k = 1:n_dim
            for l = 1:n_dim
                for p = 1:n_dim  %summation of pth eigenvector;
                    %currently using 2-D 2x2 matrix so can not do complete
                    %fix to summation of p = 1:3 later on
                    L(i,j,k,l) = L(i,j,k,l) + stor_sigma(p,1) * D(i,p)* D(j,p)* D(k,p)* D(l,p);
                    M(i,j,k,l) = M(i,j,k,l) +                   D(i,p)* D(j,p)* D(k,p)* D(l,p);
                end
            end
        end
    end
end




