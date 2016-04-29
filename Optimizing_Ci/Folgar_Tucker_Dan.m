function [DA_Dt] = Folgar_Tucker_a2_change(W,A,D, A4,CI)

gamma_dot_magnitude=sqrt(0.5*(sum(sum(gamma_dot.^2))));
Da2_Dt=zeros(3,3);


Da2_Dt= (omega*a2-a2*omega)+lambda*(gamma_dot*a2+a2*gamma_dot-2*dyadic_product(2*gamma_dot.a4)+2*C1*
end 