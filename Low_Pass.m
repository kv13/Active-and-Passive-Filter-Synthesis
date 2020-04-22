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
w_p = 2 * pi * f_p*1000
f_s = 2 * f_p
w_s = 2* pi * f_s*1000
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
n = ceil((log((10^(a_min/10)-1)/(10^(a_max/10)-1)))/ (2*log((w_s/w_p))))

%sixnotita imisias isxuos 
w_0 = w_p/((10^(a_max/10)-1)^(1/(2*n)))


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
k_f   = w_0
k_m   = 10^8 / k_f
R_1   = k_m*R_k_1



end

