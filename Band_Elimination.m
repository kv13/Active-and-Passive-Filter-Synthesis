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
bw  = w_2 - w_1;

W_p = 1;
W_s = (w_2 - w_1) / (w_4 - w_3);

W_s_k = 1;
W_p_k = 1 / W_s ; 

%taksi filtrou
n = acosh(((10^(a_min/10)-1)/(10^(a_max/10)-1))^(1/2)) / acosh(1/W_p_k);
n = ceil (n) ; 
e = 1 / sqrt(10^(a_min/10) -1);
a = asinh(1/e)/n;


%sixnotita imisias isxuos 
w_hp = 1 / cosh(1/n * acosh(1/e)); 


%poloi kai midenika sinartisis metaforas kai Q rizwn

%prokiptei n = 5 => gwnies Butterworth
y(1) = 0;
y(2) = 36; 
y(3) = 72;


for i=1:length(y)
    %upologismos polwn Chebyshev
    poles_chebyshev_real(i) = -sinh(a)*cos(y(i)*pi/180);
    poles_chebyshev_imag(i) =  cosh(a)*sin(y(i)*pi/180);
    
    W(i) = sqrt(poles_chebyshev_real(i) ^2 + poles_chebyshev_imag(i)^2);
    Q(i) = W(i)/(2*abs(poles_chebyshev_real(i)));
    
    %upologismos polwn Inverse Chebyshev
    WICH(i) = 1/W(i);
    QICH(i) = Q(i);
    
    %klimakopoiisi stin sixnotita pou mas endiaferei 
    WWICH(i) = WICH(i) * W_s;
    WQICH(i) = QICH(i)  ;
end
%midenika
Wz(1) = sec(1*pi/(2*n));
Wz(2) = sec(3*pi/(2*n));
Wz(3) = sec(5*pi/(2*n));
Wz(3)

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
    
    %antistrofi twn midenikwn
    IWWz(i) = 1 / WWz(i) ;
    
end

%GEFFE ALGORITHM 

%metasximatismos polwn
temp = 1;
qc = w_0 / bw ;
for i = 1:length(y)
   if(P(i) == 0 )
       SS(i) = abs(S(i));
       Wo(temp) = w_0;
       temp = temp + 1;
       QQ(i) = qc / SS(i);
   else
       SS(i) = abs(S(i));
       WW(i) = abs(P(i));
       
       C = SS(i)^2 + WW(i)^2;
       D = 2 * SS(i) / qc; 
       E = 4 + C / qc^2;
       G = sqrt(E^2 - 4 * D^2);
       QQ(i) = 1/D * sqrt(1/2 * ( E + G));
       K = SS(i) * QQ(i) / qc;
       W = K + sqrt ( K^2 - 1 ) ;
       Wo(temp) = W * w_0;
       temp = temp + 1 ;
       Wo(temp) = 1/W * w_0;
       temp = temp + 1;
   end
   
end

%metasximatismos midenikwn 
temp = 1;
for i=1:3
   K = 2 + IWWz(i)^2 / qc^2 ;  
   X = ( K + sqrt( K ^2 -4)) / 2;
   Wzo(temp) = w_0 * sqrt(X);
   temp = temp + 1; 
   Wzo(temp) = w_0 / sqrt(X);
   temp = temp + 1 ;
end



%%%
%%%To display use fprintf('%s will be %d \n',name , age);
%%%

end

