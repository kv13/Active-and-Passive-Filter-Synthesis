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
if(a_max < 0.1) a_max = 0.1;
W_p   = 1;
W_s   = (w_4 - w_3)/(w_2-w_1)
w_o   = sqrt(w_1*w_2);
bw    = w_2 - w_1;
q_c   = w_o / bw ;
e     = sqrt(10^(a_max/10)-1);

%ipologismos taksis filtrou
n = acosh(sqrt((10^(a_min/10) -1)/ (10^(a_max/10)- 1)))/acosh(W_s);
n = ceil(n);
a = 1/n * asinh(1/e);

%sixnotita imisias isxuos 
w_hp = cosh((1/n)*(acosh(1/e)));

%poloi chebyshev midenika kai Q

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
end

%GEFFE ALGORITHM
temp=1;
for i=1:length(y)
   Sigma = abs(poles_chebyshev_real(i));
   Omega = abs(poles_chebyshev_imag(i));
   C     = Sigma^2 +Omega^2;
   D     = 2*Sigma / q_c;
   E     = 4 + C/q_c^2;
   G     = sqrt(E^2 - 4*D^2);
   QQ    = 1/D * sqrt(1/2*(E+G));
   K     = Sigma* QQ /q_c;
   WW    = K + sqrt(K^2-1);
   w_transformed(temp) = WW* w_o;
   QQ_transformed(temp)= QQ;
   temp = temp +1;
   w_transformed(temp) = 1/WW * w_o;
   QQ_transformed(temp)= QQ;
   temp = temp +1;
   
end
fprintf('---------------------------MONADES-------------------------------\n');
fprintf('MONADA 1 Wo = %d , Q = %d \n',w_transformed(2),QQ_transformed(2));
if(QQ_transformed(2)<5)
    [numerator1,denumerator1,R1_w(1),R2_w(1),C1_w(1),C2_w(1),Z12_w(1),Z13_w(1),Hd1,RA1,RB1]=DelyiannisCircuit(QQ_transformed(2),w_transformed(2),1,w_o);
else 
    [numerator1,denumerator1,R1_w(1),R2_w(1),C1_w(1),C2_w(1),Z12_w(1),Z13_w(1),Hd1,RA1,RB1]=DelyiannisCircuit(QQ_transformed(2),w_transformed(2),2,w_o);
end
transferfunction1 = tf(numerator1,denumerator1);
transferfunction1 = (1/abs(freqresp(transferfunction1, w_o)))*transferfunction1;
plot_transfer_function(transferfunction1,[f_1 f_0 f_2 f_4 f_3] );



fprintf('MONADA 2 Wo = %d , Q = %d \n',w_transformed(3),QQ_transformed(3));
if(QQ_transformed(3)<5)
    [numerator2,denumerator2,R1_w(2),R2_w(2),C1_w(2),C2_w(2),Z12_w(2),Z13_w(2),Hd2,RA2,RB2]=DelyiannisCircuit(QQ_transformed(3),w_transformed(3),1,w_o);
else 
    [numerator2,denumerator2,R1_w(2),R2_w(2),C1_w(2),C2_w(2),Z12_w(2),Z13_w(2),Hd2,RA2,RB2]=DelyiannisCircuit(QQ_transformed(3),w_transformed(3),2,w_o);
end
transferfunction2 = tf(numerator2,denumerator2);
transferfunction2 = (1/abs(freqresp(transferfunction2, w_o)))*transferfunction2;
plot_transfer_function(transferfunction2,[f_1 f_0 f_2 f_4 f_3] );



fprintf('MONADA 3 Wo = %d , Q = %d \n',w_transformed(4),QQ_transformed(4));
if(QQ_transformed(4)<5)
    [numerator3,denumerator3, R1_w(3), R2_w(3), C1_w(3), C2_w(3), Z12_w(3), Z13_w(3),Hd3,RA3,RB3]=DelyiannisCircuit(QQ_transformed(4),w_transformed(4),1,w_o);
else
    [numerator3,denumerator3, R1_w(3), R2_w(3), C1_w(3), C2_w(3), Z12_w(3), Z13_w(3),Hd3,RA3,RB3]=DelyiannisCircuit(QQ_transformed(4),w_transformed(4),2,w_o);
end
transferfunction3 = tf(numerator3,denumerator3); 
transferfunction3 = (1/abs(freqresp(transferfunction3, w_o)))*transferfunction3;
plot_transfer_function(transferfunction3,[f_1 f_0 f_2 f_4 f_3] );



fprintf('MONADA 4 Wo = %d , Q = %d \n',w_transformed(5),QQ_transformed(5));
if(QQ_transformed(5)<5)
    [numerator4,denumerator4, R1_w(4), R2_w(4), C1_w(4), C2_w(4), Z12_w(4), Z13_w(4),Hd4,RA4,RB4]=DelyiannisCircuit(QQ_transformed(5),w_transformed(5),1,w_o);
else
    [numerator4,denumerator4, R1_w(4), R2_w(4), C1_w(4), C2_w(4), Z12_w(4), Z13_w(4),Hd4,RA4,RB4]=DelyiannisCircuit(QQ_transformed(5),w_transformed(5),2,w_o);
end
transferfunction4 = tf(numerator4,denumerator4); 
transferfunction4 = (1/abs(freqresp(transferfunction4, w_o)))*transferfunction4;
plot_transfer_function(transferfunction4,[f_1 f_0 f_2 f_4 f_3] );



fprintf('MONADA 5 Wo = %d , Q = %d \n',w_transformed(6),QQ_transformed(6));
if(QQ_transformed(6)<5)
    [numerator5,denumerator5, R1_w(5), R2_w(5), C1_w(5), C2_w(5), Z12_w(5), Z13_w(5),Hd5,RA5,RB5]=DelyiannisCircuit(QQ_transformed(6),w_transformed(6),1,w_o);
else
    [numerator5,denumerator5, R1_w(5), R2_w(5), C1_w(5), C2_w(5), Z12_w(5), Z13_w(5),Hd5,RA5,RB5]=DelyiannisCircuit(QQ_transformed(6),w_transformed(6),2,w_o);
end
transferfunction5 = tf(numerator5,denumerator5); 
transferfunction5 = (1/abs(freqresp(transferfunction5, w_o)))*transferfunction5;
plot_transfer_function(transferfunction5,[f_1 f_0 f_2 f_4 f_3] );


%Total Transfer Function
Total1 = series(transferfunction1,transferfunction2);
Total2 = series(Total1, transferfunction3); 
Total3 = series(Total2, transferfunction4);
Total  = series(Total3, transferfunction5);
plot_transfer_function(Total , [f_1 f_0 f_2 f_4 f_3]);
%Total_G =Hd1*Hd2*Hd3*Hd4*Hd5;
%Total_transfer= Total_G*Total;

%Aposvesi xwris rithmisi kerdous 
Total_Aposvesi = inv(Total);
plot_transfer_function(Total_Aposvesi , [f_1 f_0 f_2 f_4 f_3]);

%rithmisi kerdous AEM[4]=8 => Gaindb = 5dB
GaindB = 5;
Zg1 = 100;
Zg2 = (sqrt(10)-1) * Zg1;
TransferFunc = sqrt(GaindB)*Total;
plot_transfer_function(TransferFunc , [f_1 f_0 f_2 f_4 f_3]);
plot_transfer_function(inv(TransferFunc) , [f_1 f_0 f_2 f_4 f_3]);
fprintf('----TELOS----');

end

