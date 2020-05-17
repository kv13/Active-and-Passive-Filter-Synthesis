function [numerator,denumerator , R1_wanted , R2_wanted , C1_wanted , C2_wanted , Z12_wanted , Z13_wanted,Gain,RA,RB] = DelyiannisCircuit(Q,W,dec,w_o)
if(dec ==1)
    C1=1;
    C2=1;
    %Kanonikopoiisi
    Wo=1;
    R1=1/(2*Q);
    R2=2*Q;
    kf=W;
    %apo tin ekfonisi kai logo AEM(2)=5=>C=1uF
    km=10^6/W;
    C1_wanted= C1/(km*kf);
    C2_wanted=C2/(km*kf);
    R1_wanted=R1*km;
    R2_wanted=R2*km;
    T=2*Q^2;
    %apo tin ekfonisi kai logo AEM(4)=8 prepei rithmisi kerdous 5dB
    Gain_db = 5;
    Gain = 10^(Gain_db/20);
    a = 1/T;
    %rixnoume to kerdos kata a
    Z12 = R1/a;
    Z13 = R1/(1-a);
    Z12_wanted = Z12 * km;
    Z13_wanted = Z13 * km;
    RA=0;
    RB=0;
    numerator= [ -(2*W*Q) 0 ];%maybe this is wrong
    denumerator = [ 1 W/Q W^2];
elseif(dec==2)
    b = 1;
    C1= 1;
    C2= 1;
    R1= 1/sqrt(b);
    R2= sqrt(b);
    K=(Q*(b+1)-sqrt(b))/(2*Q-sqrt(b));
    kf=W;
    km=10^6/W;
    %T=K*b/(2*(K-1)-b);
    C1_wanted = C1/(km*kf);
    C2_wanted = C2/(km*kf);
    R1_wanted = R1*km;
    R2_wanted = R2*km;
    %apo tin ekfonisi kai logo AEM(4)=8 prepei rithmisi kerdous 5dB
    T=((K*w_o)/((K-1)*R1_wanted*C1_wanted))/sqrt((W^2-w_o^2)^2+(W*w_o/Q)^2);
    Gain_db = 5;
    Gain = 10^(Gain_db/20);
    a = 1/T;
    Z12= R1/a;
    Z13= R1/(1-a);
    Z12_wanted=Z12*km;
    Z13_wanted=Z13*km;
    RA=1000;
    RB=(K-1)*RA;
    numerator= [ -(2*Q*W) 0 ];%maybe this is wrong
    denumerator = [ 1 W/Q W^2];
end
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

end

