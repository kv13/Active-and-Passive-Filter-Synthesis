function [K , numerator , denumerator] = HighPass_Notch(w0,wz,Q)
%temporary values;
wt0 = w0/w0;
wtz = wz/w0;

K1 = wt0^2/wtz - 1; 
K2 = (2+K1)* Q^2 / ((2+K1)*Q^2+1);
K  = K2*(wt0^2/wtz^2);
R1 = 1;
R3 = 1;
R2 = Q^2 * (K1 + 2)^2;
R4 = Q^2 * ( K1 + 2 );
C  = 1/(Q * (2+K1));

%klimakopoiisi
k_f = w0;
%oi piknwtes apo tin ekfonisi tis ergasias exoun sygkekrimeno C=0.1μF
k_m = 10^7 * C /k_f;

%pragmatika stoixeia
R1 = k_m * R1;
R3 = R1;
R2 = R2 * k_m;
R4 = R4 * k_m;
C = 0.1;

numerator = [K 0 K*wz^2];
denumerator = [1 w0/Q w0^2];

end

