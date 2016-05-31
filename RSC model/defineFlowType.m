function [ omega, gamma_dot ] = InitializeFlowVars( flow_type,ndim )
%Initializes the vorticity and the shear rate, needed to compute the change of the orientation tensor
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

