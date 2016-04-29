function [ omega, gamma_dot ] = InitializeFlowVars(shear_rate, epsilon_dot, flow_type,ndim )

%InitializeFlowVars Initializes the vorticity and rate of strain tensor
%   in
%   flow_type   :  type of flow, 1-Simple shear XY plane vx=Gy, 2-Uniaxial
%                  elongational vx= 2Ex vy =-Ey vz = -Ez, 3-biaxial
%                  elongational vx =Ex vy=Ey vz=-2E, 4-Shearing/stretching
%                  flow vx=-Ex+Gy vy=-Ey vz=2Ez, 5-shearing flow XZ plane vx=GZ 
%   ndim        :  number of dimensions

%   out
%   gamma_dot   :  the rate of strain tensor, 3x3 matrix
%   totalStrain :  total deformation the system will undergo


nabla_v=zeros(ndim);
switch flow_type
    case (1) %simple shear
        nabla_v(1,2)=shear_rate;
    case (2) %uniaxial elongational
        nabla_v(1,1)=2*epsilon_dot;
        nabla_v(2,2)=-epsilon_dot;
        nabla_v(3,3)=-epsilon_dot;
     case (3) %biaxial elongational
        nabla_v(1,1)=epsilon_dot;
        nabla_v(2,2)=epsilon_dot;
        nabla_v(3,3)=-2*epsilon_dot;
    case (4)  %shearing/stretching flow
        nabla_v(1,2)=shear_rate;
        nabla_v(1,1)=-epsilon_dot;
        nabla_v(2,2)=-epsilon_dot;
        nabla_v(3,3)=2*epsilon_dot;
    case (5)  %shearing flow Vx = Gz
        nabla_v(1,3)=shear_rate;
end

omega=zeros(ndim);
gamma_dot=zeros(ndim);
for j=1:ndim
    for i=1:ndim
        omega(i,j)    =(nabla_v(j,i)-nabla_v(i,j));
        gamma_dot(i,j)= nabla_v(j,i)+nabla_v(i,j);
    end
end

end

