Mutconc=.25;
WTconc=.75;

load RANGE_RAND_UNSTIM;
[maxits,n]=size(RANGE_RAND_UNSTIM)

RasGTPRAND=zeros(maxits,1);
simcount=zeros(maxits,1);
differror=zeros(maxits,1);
minconc=zeros(maxits,1);
RasCons=zeros(maxits,1);
EffCons=zeros(maxits,1);

tic
parfor ii=1:maxits
    ii
    RANDfactor=RANGE_RAND_UNSTIM(ii,:);
    [a,b,count,errorout,minout,RasTotCons,EffTotCons]=ssRas_RAND(RANDfactor,Mutconc,WTconc,10);
    RasGTPRAND(ii,1)=a;
    EFFRASRAND(ii,1)=b;
    differror(ii,1)=errorout;
    simcount(ii,1)=count;
    minconc(ii,1)=minout;
    RasCons(ii,1)=RasTotCons;
    EffCons(ii,1)=EffTotCons;
end
toc

save RasGTPlevels_MillMuts_STIM