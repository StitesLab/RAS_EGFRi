load RasGTPlevels_MillMuts
load RANDfactors

delta=.05;

[m,n]=size(RasGTPRAND);
nomuts=max(m,n);

[CommonMuts_RasGTP(1),CommonMuts_EffRasGTP(1)]=ssRas(1,0.25,0.75,1);
[CommonMuts_RasGTP(2),CommonMuts_EffRasGTP(2)]=ssRas(2,0.25,0.75,1);
[CommonMuts_RasGTP(3),CommonMuts_EffRasGTP(3)]=ssRas(3,0.25,0.75,1);

min_RASGTP=min(CommonMuts_RasGTP)-delta;
max_RASGTP=max(CommonMuts_RasGTP)+delta;

muts_in_range=0;

aaa=find(RasGTPRAND>min_RASGTP);
RAS_2=RasGTPRAND(aaa);
EFF_2=EFFRASRAND(aaa);
RANDfactors_2=RANDfactors(aaa,:);

bbb=find(RAS_2<max_RASGTP);
RAS_3=RAS_2(bbb);
EFF_3=EFF_2(bbb);
RANDfactors_3=RANDfactors_2(bbb,:);

min_EffRas=min(CommonMuts_EffRasGTP)-delta;
max_EffRas=max(CommonMuts_EffRasGTP)+delta;

ccc=find(EFF_3>min_EffRas);
RAS_4=RAS_3(ccc);
EFF_4=EFF_3(ccc);
RANDfactors_4=RANDfactors_3(ccc,:);

ddd=find(EFF_4<max_EffRas);
RAS_IN_RANGE_UNSTIM=RAS_4(ddd);
EFF_IN_RANGE_UNSTIM=EFF_4(ddd);
RANGE_RAND_UNSTIM=RANDfactors_4(ddd,:);


save RANGE_RAND_UNSTIM RANGE_RAND_UNSTIM
save RAS_IN_RANGE_UNSTIM RAS_IN_RANGE_UNSTIM
save EFF_IN_RANGE_UNSTIM EFF_IN_RANGE_UNSTIM