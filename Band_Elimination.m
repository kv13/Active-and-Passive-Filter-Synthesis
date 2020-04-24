function [] = Band_Elimination()

AEM = [8 5 0 8];

%prodiagrafes
f_0 = 1800;
f_1 = 1200 +25 * (9 - AEM(4));
f_2 = f_0^2 /f_1;
D   = (1/1.8) * (f_0^2 - f_1^2)/f_1;

f_3 = (-D + sqrt(D^2 + 4*f_0^2))/2;
f_4 = f_0^2 / f_3;
a_min = 30 - AEM(3);
a_max = 0.5 + AEM(4)/18;


w_0 = 2 * pi * f_0;
w_1 = 2 * pi * f_1;
w_2 = 2 * pi * f_2;
w_3 = 2 * pi * f_3;
w_4 = 2 * pi * f_4;


W_p = 1;
W_s = (w_2 - w_1) / (w_4 - w_3);

W_s_k = 1;
W_p_k = 1 / W_s ; 

%taksi filtrou
n = acosh(((10^(a_min/10)-1)/(10^(a_max/10)-1))^(1/2)) / acosh(1/W_p_k);
n = ceil (n) ; 
e = 1 / sqrt(10^2.5 -1);
a = 1/n * asinh(1/e);


%sixnotita imisias isxuos 
w_hp = 1 / cosh(1/n * acosh(1/e)); 


%poloi kai midenika sinartisis metaforas kai Q rizwn

%prokiptei n = 5 => gwnies Butterworth
y(1) = 0;
y(2) = 36; 
y(3) = 72;


for i=1:length(y)
    %upologismos polwn Chebyshev
    poles_chebyshev_real(i) = -sinh(a)*cos(y(i));
    poles_chebyshev_imag(i) =  cosh(a)*sin(y(i));
    
    W(i) = sqrt(poles_chebyshev_real(i) ^2 + poles_chebyshev_imag(i)^2);
    Q(i) = W(i)/(2*abs(poles_chebyshev_real(i)));
    
    %upologismos polwn Inverse Chebyshev
    WICH(i) = 1/W(i);
    QICH = Q(i);
    
    %klimakopoiisi stin sixnotita pou mas endiaferei 
    WWICH(i) = WICH(i) * W_s;
    WQICH(I) = QICH(i)  ;
end
%midenika
Wz(1) = sec(1*pi/(2*n));
Wz(2) = sec(3*pi/(2*n));
Wz(3) = sec(5*pi/(2*n));

%klimakopoiisi midenikwn
WWz(1) = Wz(1) * W_s;
WWz(2) = Wz(2) * W_s;
WWz(3) = Wz(3) * W_s;

%poloi anwdiavatis sinartisis 
for i = 1:length(y)
    %antistrofi polwn ICH
    X(i) = 1 / WWICH(i);
    S(i) =- X(i) / (2* WQICH(i));
    P(i) = sqrt(X(i)^2 - S(i)^2);
    
end




end

