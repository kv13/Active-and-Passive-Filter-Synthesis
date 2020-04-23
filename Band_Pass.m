%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%               KONSTANTINOS VERGOPOULOS               %%%%%%%
%%%%%%%          AEM 8508 MAIL:vkonstant@ece.auth.gr         %%%%%%%
%%%%%%%             BAND PASS FILTER : CHEBYSHEV             %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = Band_Pass()

AEM = [8 5 0 8];

%prodiagrafes
f_0 = 0.65 * 1000;
f_1 = 400 + 25 * AEM(3);
f_2 = f_0^2 / f_1 ;
D   = 2.3 * (f_0^2 - f_1^2)/f_1;
f_3 = ( -D + sqrt(D^2 + 4*f_0^2))/2;
f_4 = f_0^2/f_3;

w_0 = 2 * pi * f_0;
w_1 = 2 * pi * f_1;
w_2 = 2 * pi * f_2;
w_3 = 2 * pi * f_3;
w_4 = 2 * pi * f_4;

a_min = 27.5 + AEM(4);
a_max = 0.5 + (AEM(3)-5)/10;
W_p   = 1;
W_s   = (w_4 - w_3)/(w_2-w_1)
w_o   = sqrt(w_1*w_2);
bw    = w_2 - w_1;
e     = sqrt(10^(a_max/10)-1);

%ipologismos taksis filtrou
n = acosh(sqrt((10^(a_min/10) -1)/ (10^(a_max/10)- 1)))/acosh(W_s)
a = 1/n * asinh(1/e);

%sixnotita imisias isxuos 
w_hp = cosh(1/(n * sqrt(acosh(10^(a_max/10) - 1))))

%poloi midenika kai Q 
end

