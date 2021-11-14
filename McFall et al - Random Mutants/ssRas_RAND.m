function [fractRasact,fractEffbound,fractWTRasact,fractMutRasact,fractEffboundWT,fractEffboundMut,count]=ssRas_RAND(RANDfactor,Mutconc,WTconc,GEFfact);

GTot=4e-7;
EffTot=4e-7;
GEF=2e-10;
GAP=6e-11;
GEF=GEF*GEFfact;

GTP=180e-6;
GDP=18e-6;
Kd=80e-9;
volscale=250;

GTPaseModuleParameters;

k(1)=GEF*kD; %VmaxD=k(1);
k(2)=GEF*kT; %VmaxT=k(2);
k(3)=KmD; %=k(3);
k(4)=KmT; %=k(4);
k(5)=GAP*kcat; %Vmax=k(5);
k(6)=Km; %=k(6);
k(7)=kint; %=k(7);
k(8)=kdissD;% =k(8);
k(9)=kdissT;%=k(9);
k(10)=kassD; %k(10);
k(11)=kassT; %k(11);
k(12)=kassEff;%=k(12);
k(13)=kdissEff;%=k(13);

CancerRANDspec;
       
%Parameters for the Ras mutant
k(14)=GEF*mkD;
k(15)=GEF*mkT;
k(16)=mKmD;
k(17)=mKmT; %V=k(17);
k(18)=GAP*mkcat; %VmaxV=k(18);
k(19)=mKm; %=k(19);
k(20)=mkint; %kintV=k(20);
k(21)=mkdissD; %kdissDV=k(21);
k(22)=mkdissT; %kdissTV=k(22);
k(23)=mkassD; %kassDGDPV=k(23);
k(24)=mkassT; %kassTGTPV=k(24);
k(25)=mkassEff; %kassEffV=k(25);
k(26)=mkdissEff;

%Initial Conditions
MutTot=Mutconc*GTot;
WTRasTot=WTconc*GTot;

y(1)=WTRasTot; %GD=y(1);
y(2)=0; %GT=y(2);
y(3)=0; %G0=y(3);
y(4)=EffTot; %Eff=y(4);
y(5)=0; %GTEff=y(5);
y(6)=MutTot; %GDV=y(6);
y(7)=0; %GTV=y(7);
y(8)=0; %G0V=y(8);
y(9)=0; %GTVEff=y(9);

yprev=zeros(1,9);
yf=y;

tspan=[0 10000];
options=odeset('RelTol',1e-9,'AbsTol',1e-9,'NonNegative',[1:9]);
count=0;
dmet=ones(1,13);

while dot(dmet,dmet)>1e-40
    yprev=yf;
    [t,y]=ode15s(@fullRasmodelwithtwoRasPools,tspan,yf,options,k);
    [m,n]=size(y);
    yf=y(m,:);
    %count=count+1
    %figure(count);
    %plot(t,y);
    dmet=(yf-yprev);
    f0difference=dot(dmet,dmet);
end


   
[m,n]=size(y);

freeGT=y(m,2);
GTEffcomplex=y(m,5);
freeGTV=y(m,7);
GTVEffcomplex=y(m,9);

fractRasact=(freeGT+freeGTV+GTEffcomplex+GTVEffcomplex)/(WTRasTot+MutTot);
fractEffbound=(GTEffcomplex+GTVEffcomplex)/EffTot;
fractWTRasact=(freeGT+GTEffcomplex)/(WTRasTot);
fractMutRasact=(freeGTV+GTVEffcomplex)/(MutTot);
fractEffboundWT=(GTEffcomplex)/EffTot;
fractEffboundMut=(GTVEffcomplex)/EffTot;