%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%               KONSTANTINOS VERGOPOULOS               %%%%%%%        
%%%%%%%          AEM 8508 MAIL:vkonstant@ece.auth.gr         %%%%%%%
%%%%%%%             LOW PASS FILTER : BUTTERWORTH            %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = Low_Pass()

%SPECIFICATIONS
AEM = [8 5 0 8];
m = 3;
f_p = 0.6 * (3+m)
w_p = 2 * pi * f_p*1000;
f_s = 2 * f_p;
w_s = 2* pi * f_s*1000;
if(AEM(4)>1)
    a_min = 17.5 + (AEM(4)-5) * 0.5;
else
    a_min = 17.5 + (1-5) * 0.5;
end
if(AEM(3)>1)
    a_max = 0.6 +(AEM(3)-3)/10;
else
    a_max = 0.6 + (1-3)/10;
end

%taksi filtorou 
n = ceil((log((10^(a_min/10)-1)/(10^(a_max/10)-1)))/ (2*log((w_s/w_p))));

%sixnotita imisias isxuos 
w_0 = w_p/((10^(a_max/10)-1)^(1/(2*n)));


y(1) = 0;
y(2) = 36;
y(3) = 72;

for i=1:3
    realPart(i) = cos(y(i)*(pi/180));
    imPart(i)   = sin(y(i)*(pi/180));
    Q(i)        = 1/(2*cos(y(i)*(pi/180)));
end

%1h monada
R_k_1 = 1 ;
C_k_1 = 1;
%1h monada klimakopoiisi
k_f   = w_0;
k_m   = 10^8 / k_f;
R_1   = k_m*R_k_1;
c_1   = 1/(k_m*k_f);


%2h monada Sallen-Key strategy 1
R_k_21 = 1 ;
R_k_22 = 1 ;
C_k_21 = 1 ;
C_k_22 = 1 ;
w_20   = 1 ;
k_2    = 3 - 1/Q(2);
r_k_21   = 1 ;
r_k_22   = 2 - 1/Q(2);

%2h monada klimakopoiish
C_21 = C_k_21/(k_f * k_m);
C_22 = C_k_22/(k_f *k_m);
R_21 = R_k_21 * k_m;
R_22 = R_k_22 * k_m;
r_21 = r_k_21 * k_m;
r_22 = r_k_22 * k_m;


%3h monada Sallen-Key strategy 1
R_k_31 = 1 ;
R_k_32 = 1 ;
C_k_31 = 1 ;
C_k_32 = 1 ;
w_30   = 1 ;
k_3    = 3 - 1/Q(3);
r_k_31   = 1 ;
r_k_32   = 2 - 1/Q(3);

%3h monada klimakopoiish
C_31 = C_k_31/(k_f * k_m);
C_32 = C_k_32/(k_f *k_m);
R_31 = R_k_31 * k_m;
R_32 = R_k_32 * k_m;
r_31 = r_k_31 * k_m;
r_32 = r_k_32 * k_m;

%Rithimisi Kerdous
Gain_db = 0 ;
Gain = 10^( Gain_db / 20 );
Gain_total = 1 * k_2 * k_3 ; 
a = Gain / Gain_total ;
Z_2 = R_21 / a;
Z_3 = R_21 / (1-a);

% sinartisis metaforas
f_p_hz = f_p *1000;
f_s_hz = f_s *1000;
f_0_hz = w_0 /(2*pi);
T1 = tf([w_0],[1 w_0]);
T2 = k_2*tf(w_0^2,[1 w_0/Q(2) w_0^2]);
T3 = k_3*tf(w_0^2,[1 w_0/Q(3) w_0^2]);
plot_transfer_function(T1, [f_s_hz f_0_hz f_p_hz]);
plot_transfer_function(T2, [f_s_hz f_0_hz f_p_hz]);
plot_transfer_function(T3, [f_s_hz f_0_hz f_p_hz]);

T12 = series(T1 , T2);
Tol = series(T12 , T3);
T = Tol * a;
plot_transfer_function(Tol, [f_s_hz f_0_hz f_p_hz]);
plot_transfer_function(T, [f_s_hz f_0_hz f_p_hz]);


invH = inv(Tol);
invHk = inv(T);
plot_transfer_function(invH , [f_s_hz f_0_hz f_p_hz]);
plot_transfer_function(invHk , [f_s_hz f_0_hz f_p_hz]);


end

