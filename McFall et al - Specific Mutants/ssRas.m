function [fractRasact,fractEffbound,fractWTRasact,fractMutRasact]=ssRas(Mutflag,Mutconc,WTconc,GEFfact);

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

if Mutflag==0;
    CancerG12V; 
elseif Mutflag==1;
    CancerG12V;
elseif Mutflag==2;
    CancerG12D;
elseif Mutflag==3;
    CancerG13D;

elseif Mutflag==4;
    CancerQ61H;
elseif Mutflag==5;
    CancerQ61K;
elseif Mutflag==6;
    CancerQ61L;
elseif Mutflag==7;
    CancerQ61P;
elseif Mutflag==8;
    CancerQ61R;
elseif Mutflag==9;
    CancerQ61W;

elseif Mutflag==101;
    CancerQ61K_WT;
elseif Mutflag==102;
    CancerQ61R_WT;

elseif Mutflag==444;
    CancerG13DRabara;

else
    display('Ras Mutant Undefined');
end
       
%mutant RAS parameters
k(14)=GEF*mkD;
k(15)=GEF*mkT;
k(16)=mKmD;
k(17)=mKmT;
k(18)=GAP*mkcat;
k(19)=mKm;
k(20)=mkint;
k(21)=mkdissD;
k(22)=mkdissT;
k(23)=mkassD;
k(24)=mkassT;
k(25)=mkassEff;
k(26)=mkdissEff; 

MutTot=Mutconc*GTot;
WTRasTot=WTconc*GTot;

y(1)=WTRasTot;
y(2)=0;
y(3)=0; 
y(4)=EffTot;
y(5)=0; 
y(6)=MutTot;
y(7)=0; 
y(8)=0; 
y(9)=0; 

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