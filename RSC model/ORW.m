function [ a4 ] = ORW( a2 )
% Computes the Orthotropic  Fittedclosure approximation
% the coefficients are from Chung and Kwon (2001) 
a4 =zeros(3,3,3,3);

C = [0.070055  0.339376 0.590331 -0.396796  0.333693 0.411944; ...
     0.115177 -0.368267 0.252880  0.094820  0.800181 0.535224; ...
     1.249811 -2.148297 0.898521 -2.290157  1.044147 1.934914];

% C = [0.060964 0.371243 0.555301 -0.369160 0.318266 0.371218;    ...
%     0.124711 -0.389402 0.258844 0.086169 0.796080 0.544992;   ...
%     1.228982 -2.054116 0.821548 -2.260574 1.053907 1.819756];

%Compute A11 A22 and A33 of the othotropic tensor
A = C * [1 ; a2(1,1) ; a2(1,1)*a2(1,1) ; a2(2,2) ; a2(2,2)*a2(2,2) ; a2(1,1)*a2(2,2) ];

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

