function [R1 , R2 , R3 , R4 , R5 , C, K] = LowPass_Notch(w0 , wz , Q)
%temporary values
wt0 = w0/w0;
wtz = wz/w0;

R1 = 1;
R4 = 1;
C  = 1/(2*Q);
R2 = 4* Q^2;
R5 = 4*Q^2 / (wtz^2 -1);
K  = 1/(1+ wtz^2/(2*Q^2));
R3 = wtz^2/(2*Q^2);

%klimakopoiisi 
k_f = w0;
%oi piknwtes apo tin ekfonisi tis ergasias exoun sygkekrimeno C=0.1μF
k_m = 10^7 * C /k_f;
R1 = R1 * k_m;
R2 = R2 * k_m;
R3 = R3 * k_m;
R4 = R4 * k_m;
R5 = R5 * k_m;
C= 0.1;



end

