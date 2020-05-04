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

%%ektipwsi sixnotitwn f
fprintf('-----------------------------------------------------------------\n');
fprintf('f0 = %d \n',f_0);
fprintf('f1 = %d \n',f_1);
fprintf('f2 = %d \n',f_2);
fprintf('f3 = %d \n',f_3);
fprintf('f4 = %d \n',f_4);
fprintf('-----------------------------------------------------------------\n');


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
%ektipwsi kiklikwn sixnotitwn
fprintf('-----------------------------------------------------------------\n');
fprintf('w0 = %d \n',w_0);
fprintf('w1 = %d \n',w_1);
fprintf('w2 = %d \n',w_2);
fprintf('w3 = %d \n',w_3);
fprintf('w4 = %d \n',w_4);
fprintf('Wp = %d \n',W_p);
fprintf('Ws = %d \n',W_s);
fprintf('kanonikopoiisi :\n');
fprintf('W_p_k = %d \n',W_p_k);
fprintf('W_s_k = %d \n',W_s_k);
fprintf('-----------------------------------------------------------------\n');

%taksi filtrou
n = acosh(((10^(a_min/10)-1)/(10^(a_max/10)-1))^(1/2)) / acosh(1/W_p_k);
n = ceil (n) ; 
e = 1 / sqrt(10^(a_min/10) -1);
a = asinh(1/e)/n;


%sixnotita imisias isxuos 
w_hp = 1 / cosh(1/n * acosh(1/e)); 


%ektiposi upoloipwn stoixeiwn
fprintf('-----------------------------------------------------------------\n');
fprintf('D = %d \n',D);
fprintf('amin = %d \n',a_min);
fprintf('amax = %d \n',a_max);
fprintf('bw = %d \n',bw);
fprintf('w4 = %d \n',w_4);
fprintf('n = %d \n',n);
fprintf('e = %d \n',e);
fprintf('a = %d \n',a);
fprintf('whp = %d \n',w_hp);
fprintf('-----------------------------------------------------------------\n');

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
fprintf('-----------------------CHEBYSHEV POLES---------------------------\n');
for i=1:3
    fprintf('%d +- i%d ===> W = %d , Q = %d \n ',poles_chebyshev_real(i) , poles_chebyshev_imag(i),W(i),Q(i));
    
end
fprintf('-----------------------------------------------------------------\n');
fprintf('\n \n');
fprintf('-----------------------In.CHEBYSHEV POLES------------------------\n');
fprintf('InW = 1/W , InQ = Q\n');
for i=1:3
    fprintf('xwris klimakopoiisi: W = %d , Q = %d \n',WICH(i),QICH(i));
    fprintf('me klimakopoiisi Ws: W = %d , Q = %d \n',WWICH(i),WQICH(i));
end
fprintf('-----------------------------------------------------------------\n');
fprintf('--------------------------MIDENIKA-------------------------------\n');
for i=1:3
    fprintf('xwris klimakopoiisi: Wz = %d  \n',Wz(i));
    fprintf('me klimakopoiisi Ws: Wz = %d  \n',WWz(i));
end
fprintf('-----------------------------------------------------------------\n');

fprintf('----------------------Metatropi LP =>HP--------------------------\n');
for i=1:3
    fprintf('polos: %d +- i%d , Q=%d \n',S(i),P(i),WQICH(i));
    fprintf('mideniko: Wz = %d\n',IWWz(i));
end
fprintf('-----------------------------------------------------------------\n');
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
fprintf('------------------------BAND ELIMINATION-------------------------\n');
temp = 1;
for i=1:3
    fprintf('polos: Wo = %d , Q =%d \n',Wo(temp),QQ(i));
    fprintf('mideniko Wz = %d \n', Wzo(6-temp));
    temp = temp +1;
    if (i>1)
        fprintf('polos: Wo = %d , Q =%d \n',Wo(temp),QQ(i));
        fprintf('mideniko Wz = %d \n', Wzo(6-temp));
        temp = temp+1;
    end
    
end
fprintf('-----------------------------------------------------------------\n');
fprintf('---------------------------MONADES-------------------------------\n');
%gia tis monades isxuei oti an Wo>=Wz => HPN an Wo<Wz => LPN
fprintf('MONADA 1 Wo = %d , Q = %d , Wz = %d \n',Wo(1),QQ(1),Wzo(5));
[GAIN(1) , num1 , den1 ] = HighPass_Notch(Wo(1), Wzo(5), QQ(1));
Transferfunction1 = tf(num1 , den1);
plot_transfer_function( Transferfunction1, [f_1 f_3 f_0 f_4 f_2] );

fprintf('MONADA 2 Wo = %d , Q = %d , Wz = %d \n',Wo(2),QQ(2),Wzo(3));
[GAIN(2) , num2 , den2] = HighPass_Notch(Wo(2) , Wzo(3) , QQ(2));
Transferfunction2 = tf(num2 , den2);
plot_transfer_function( Transferfunction2, [f_1 f_3 f_0 f_4 f_2] );

fprintf('MONADA 3 Wo = %d , Q = %d , Wz = %d \n',Wo(3),QQ(2),Wzo(4));
[GAIN(3) , num3 , den3] = LowPass_Notch(Wo(3) , Wzo(4) , QQ(2));
Transferfunction3 = tf(num3 , den3);
plot_transfer_function( Transferfunction3, [f_1 f_3 f_0 f_4 f_2] );

fprintf('MONADA 4 Wo = %d , Q = %d , Wz = %d \n',Wo(4),QQ(3),Wzo(1));
[GAIN(4) , num4 , den4] = HighPass_Notch( Wo(4) , Wzo(1) , QQ(3));
Transferfunction4 = tf(num4 , den4);
plot_transfer_function( Transferfunction4, [f_1 f_3 f_0 f_4 f_2] );

fprintf('MONADA 5 Wo = %d , Q = %d , Wz = %d \n',Wo(5),QQ(3),Wzo(2));
[GAIN(5) , num5 , den5] = LowPass_Notch(Wo(5) , Wzo(2) , QQ(3));
Transferfunction5 = tf(num5 , den5);
plot_transfer_function( Transferfunction5, [f_1 f_3 f_0 f_4 f_2] );

%Total Transfer Function
Total1 = series(Transferfunction1,Transferfunction2);
Total2 = series(Total1, Transferfunction3); 
Total3 = series(Total2, Transferfunction4);
Total  = series(Total3, Transferfunction5);
plot_transfer_function(Total , [f_1 f_3 f_0 f_4 f_2] );
end

