load RasGTPlevels_MillMuts_STIM
load RANGE_RAND_UNSTIM
load RAS_IN_RANGE_UNSTIM
load EFF_IN_RANGE_UNSTIM

delta=.05;
deltaspan=.05;

[m,n]=size(RasGTPRAND);
nomuts=max(m,n);

[CommonMuts_RasGTP(1),CommonMuts_EffRasGTP(1)]=ssRas(1,0.25,0.75,10);
[CommonMuts_RasGTP(2),CommonMuts_EffRasGTP(2)]=ssRas(2,0.25,0.75,10);
[CommonMuts_RasGTP(3),CommonMuts_EffRasGTP(3)]=ssRas(3,0.25,0.75,10);

min_RASGTP=min(CommonMuts_RasGTP)-delta;
max_RASGTP=max(CommonMuts_RasGTP)+delta;

muts_in_range=0;

aaa=find(RasGTPRAND>min_RASGTP);
RAS_2=RasGTPRAND(aaa);
EFF_2=EFFRASRAND(aaa);
RAS_IN_RANGE_UNSTIM_2=RAS_IN_RANGE_UNSTIM(aaa);
EFF_IN_RANGE_UNSTIM_2=EFF_IN_RANGE_UNSTIM(aaa);
RANDfactors_2=RANGE_RAND_UNSTIM(aaa,:);

bbb=find(RAS_2<max_RASGTP);
RAS_3=RAS_2(bbb);
EFF_3=EFF_2(bbb);
RAS_IN_RANGE_UNSTIM_3=RAS_IN_RANGE_UNSTIM_2(bbb);
EFF_IN_RANGE_UNSTIM_3=EFF_IN_RANGE_UNSTIM_2(bbb);
RANDfactors_3=RANDfactors_2(bbb,:);

min_EffRas=min(CommonMuts_EffRasGTP)-delta;
max_EffRas=max(CommonMuts_EffRasGTP)+delta;

ccc=find(EFF_3>min_EffRas);
RAS_4=RAS_3(ccc);
EFF_4=EFF_3(ccc);
RAS_IN_RANGE_UNSTIM_4=RAS_IN_RANGE_UNSTIM_3(ccc);
EFF_IN_RANGE_UNSTIM_4=EFF_IN_RANGE_UNSTIM_3(ccc);
RANDfactors_4=RANDfactors_3(ccc,:);

ddd=find(EFF_4<max_EffRas);
RAS_IN_RANGE_5=RAS_4(ddd);
EFF_IN_RANGE_5=EFF_4(ddd);
RAS_IN_RANGE_UNSTIM_5=RAS_IN_RANGE_UNSTIM_4(ddd);
EFF_IN_RANGE_UNSTIM_5=EFF_IN_RANGE_UNSTIM_4(ddd);
RANDfactors_5=RANDfactors_4(ddd,:);

[U_CommonMuts_RasGTP(1),U_CommonMuts_EffRasGTP(1)]=ssRas(1,0.25,0.75,1);
[U_CommonMuts_RasGTP(2),U_CommonMuts_EffRasGTP(2)]=ssRas(2,0.25,0.75,1);
[U_CommonMuts_RasGTP(3),U_CommonMuts_EffRasGTP(3)]=ssRas(3,0.25,0.75,1);

Span_CommonMuts_RasGTP=(CommonMuts_RasGTP-U_CommonMuts_RasGTP);
Span_CommonMuts_EffRasGTP=(CommonMuts_EffRasGTP-U_CommonMuts_EffRasGTP);

minspanRas=min(Span_CommonMuts_RasGTP)-deltaspan;
maxspanRas=max(Span_CommonMuts_RasGTP)+deltaspan;
minspanEff=min(Span_CommonMuts_EffRasGTP)-deltaspan;
maxspanEff=max(Span_CommonMuts_EffRasGTP)+deltaspan;

SPAN_MUTS_RASGTP=RAS_IN_RANGE_5 - RAS_IN_RANGE_UNSTIM_5;
SPAN_MUTS_EFFRASGTP=EFF_IN_RANGE_5 - EFF_IN_RANGE_UNSTIM_5;

eee=find(SPAN_MUTS_RASGTP<maxspanRas);
SPAN_RAS_2=SPAN_MUTS_RASGTP(eee);
SPAN_EFF_2=SPAN_MUTS_EFFRASGTP(eee);
RANDfactors_6=RANDfactors_5(eee,:);

fff=find(SPAN_RAS_2>minspanRas);
SPAN_RAS_3=SPAN_RAS_2(fff);
SPAN_EFF_3=SPAN_EFF_2(fff);
RANDfactors_7=RANDfactors_6(fff,:);

ggg=find(SPAN_EFF_3<maxspanEff);
SPAN_RAS_4=SPAN_RAS_3(ggg);
SPAN_EFF_4=SPAN_EFF_3(ggg);
RANDfactors_8=RANDfactors_7(ggg,:);

hhh=find(SPAN_EFF_4>minspanEff);
SPAN_RAS_5=SPAN_RAS_4(hhh);
SPAN_EFF_5=SPAN_EFF_4(hhh);
RANGE_RAND=RANDfactors_8(hhh,:);

save RANGE_RAND RANGE_RAND
    