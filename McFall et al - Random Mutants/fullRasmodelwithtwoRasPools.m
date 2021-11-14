function dydt=fullRasmodelwithtwoRasPools(t,y,k);

GD=y(1);
GT=y(2);
G0=y(3);
Eff=y(4);
GTEff=y(5);
GDV=y(6);
GTV=y(7);
G0V=y(8);
GTEffV=y(9);


VmaxD=k(1);
VmaxT=k(2);
KmD=k(3);
KmT=k(4);
Vmax=k(5);
Km=k(6);
kint=k(7);
kdissD=k(8);
kdissT=k(9);
kassDGDP=k(10);
kassTGTP=k(11);
kassEff=k(12);
kdissEff=k(13);

VmaxDV=k(14);
VmaxTV=k(15);
KmDV=k(16);
KmTV=k(17);
VmaxV=k(18);
KmV=k(19);
kintV=k(20);
kdissDV=k(21);
kdissTV=k(22);
kassDGDPV=k(23);
kassTGTPV=k(24);
kassEffV=k(25);
kdissEffV=k(26);

R1=(VmaxD*GD/KmD-VmaxT*GT/KmT)/(1+GD/KmD+GT/KmT+GDV/KmDV+GTV/KmTV);
R2=Vmax*GT/(Km*(1+GTV/KmV)+GT);
R3=kint*GT;
R4=kdissD*GD-kassDGDP*G0;
R5=kdissT*GT-kassTGTP*G0;
R6=kassEff*GT*Eff-kdissEff*GTEff;
R7=kint*GTEff;

R8=(VmaxDV*GDV/KmDV-VmaxTV*GTV/KmTV)/(1+GD/KmD+GT/KmT+GDV/KmDV+GTV/KmTV);
R9=VmaxV*GTV/(KmV*(1+GT/Km)+GTV);
R10=kintV*GTV;
R11=kdissDV*GDV-kassDGDPV*G0V;
R12=kdissTV*GTV-kassTGTPV*G0V;
R13=kassEffV*GTV*Eff-kdissEffV*GTEffV;
R14=kintV*GTEffV;


dydt=[-R1+R2+R3-R4+R7; %GD
    R1-R2-R3-R5-R6; %GT
    R4+R5; %G0
    -R6+R7-R13+R14; %Eff
    R6-R7; %GTEff
     -R8+R9+R10-R11+R14;% GDV
     R8-R9-R10-R12-R13;%GTV
     R11+R12;%G0V
     R13-R14];%GTEffV