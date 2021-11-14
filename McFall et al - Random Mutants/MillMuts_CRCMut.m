Mutconc=.25; 
WTconc=.75;

load RANDfactors;

[maxits,n]=size(RANDfactors)

RasGTPRAND=zeros(maxits,1);
simcount=zeros(maxits,1);
differror=zeros(maxits,1);
minconc=zeros(maxits,1);
RasCons=zeros(maxits,1);
EffCons=zeros(maxits,1);

parfor ii=1:maxits
    ii
    RANDfactor=RANDfactors(ii,:);
    [a,b,count,errorout,minout,RasTotCons,EffTotCons]=ssRas_RAND(RANDfactor,Mutconc,WTconc,1);
    RasGTPRAND(ii,1)=a;
    EFFRASRAND(ii,1)=b;
    differror(ii,1)=errorout;
    simcount(ii,1)=count;
    minconc(ii,1)=minout;
    RasCons(ii,1)=RasTotCons;
    EffCons(ii,1)=EffTotCons;
end
save RasGTPlevels_MillMuts