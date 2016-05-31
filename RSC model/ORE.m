function [ a4 ] = ORE( a2 )
% Computes the Orthotropic  Fittedclosure approximation

% find eginevalues of a2
%lambda = eig(a2);

% fitted values by Cintra
a4 =zeros(3,3,3,3);

C = [  0.636256797 -1.872662964 -4.479708732 11.95895623 3.844596924 11.34209243 -10.95826261 -20.72779947 ...
        -2.116232145 -12.38756329 9.815983897 3.479015106 11.74929112 0.508041387 4.883665978 ;...
        0.636256797 -3.315272297 -3.037099398 11.8273286 6.881539521 8.436777 -15.91206672 -15.15158726 ...
        -6.487289336 -8.638914193 9.325203435 7.746837517 7.481468706 2.284765316 3.597722511 ; ...
        2.740532896 -9.121965098 -12.2570587 34.31990189 13.82946991 25.86847553 -37.7029118 -50.27564319 ...
        -10.88017611 -26.96369152 27.33467981 15.26506861 26.1134914 3.432138403 10.61174181];
lambda= [a2(1,1) a2(2,2) a2(3,3)  ] ;

vect =  [1; lambda(1); lambda(2); lambda(1)*lambda(2); lambda(1)^2;lambda(2)^2;...
          lambda(1)^2*lambda(2); lambda(2)^2*lambda(1);lambda(1)^3;lambda(2)^3; lambda(1)^2*lambda(2)^2;...
         lambda(1)^3*lambda(2); lambda(2)^3*lambda(1); lambda(1)^4 ; lambda(2)^4   ] ;
 

%Compute A11 A22 and A33 of the othotropic tensor
%A = C * [1 ; a2(1,1) ; a2(1,1)*a2(1,1) ; a2(2,2) ; a2(2,2)*a2(2,2) ; a2(1,1)*a2(2,2) ];
A = C * vect;
A2 = [0 1 1; 1 0 1; 1 1 0 ] \ ( [a2(1,1);a2(2,2);a2(3,3)  ] - A  );

A_mn = [  A(1) A2(3) A2(2)   0      0    0; ...
    A2(3)  A(2) A2(1)   0      0    0 ;...
    A2(2) A2(1)  A(3)   0      0    0 ;...
    0     0     0    A2(1)    0    0 ;...
    0     0     0      0    A2(2)  0; ...
    0     0     0      0      0   A2(3)];

for i = 1 :3
    for j = 1 :3
        switch i+j
            case 2
                m = 1;
            case 4
                if i == 1 || i== 3
                    m = 5;
                else
                    m = 2;
                end
            case 6
                m = 3;
            case 5
                m = 4;
            case 3
                m = 6;
        end
        for k = 1 :3
            for l =1 :3
                switch k+l
                    case 2
                        n = 1;
                    case 4
                        if k == 1 || k== 3
                            n = 5;
                        else
                            n = 2;
                        end
                    case 6
                        n = 3;
                    case 5
                        n = 4;
                    case 3
                        n = 6;
                end              
                a4(i,j,k,l) = A_mn(m,n);     
            end
        end
    end
end



end

