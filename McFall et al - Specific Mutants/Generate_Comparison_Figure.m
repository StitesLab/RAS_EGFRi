GEFvals=[10:-.25:1];
[m,n]=size(GEFvals);

RasIndices=[0,1,2,3,444];

for ii=1:n
    GEFfact=GEFvals(ii);

    for jj=1:5
    
        mutflag=RasIndices(jj);
        if jj==1;
            [RasGTP(ii,jj),RasEff(ii,jj),WTRAS(ii,jj),MUTRAS(ii,jj)]=ssRas(mutflag,0,1,GEFfact);
        else
            [RasGTP(ii,jj),RasEff(ii,jj),WTRAS(ii,jj),MUTRAS(ii,jj)]=ssRas(mutflag,.25,.75,GEFfact);
        end
    end
end

RasGTP=RasGTP*100;
RasEff=RasEff*100;


figure(751)
hold on
plot(GEFvals,100*RasGTP(:,1)/(RasGTP(1,1)),'k');
plot(GEFvals,100*RasGTP(:,2)/(RasGTP(1,2)),'g');
plot(GEFvals,100*RasGTP(:,3)/(RasGTP(1,3)),'b');
plot(GEFvals,100*RasGTP(:,4)/(RasGTP(1,4)),'r'); %G13D - Mech 1
plot(GEFvals,100*RasGTP(:,5)/(RasGTP(1,5)),'m'); %G13D - Mech 2
set(gca, 'XDir','reverse');
set(gca,'XScale','log');
xlabel('(EGFR) SOS1/2 Activity');
ylabel('RAS-GTP (% max)');

figure(752)
hold on
plot(GEFvals,100*WTRAS(:,1)/(WTRAS(1,1)),'k');
plot(GEFvals,100*WTRAS(:,2)/(WTRAS(1,2)),'g');
plot(GEFvals,100*WTRAS(:,3)/(WTRAS(1,3)),'b');
plot(GEFvals,100*WTRAS(:,4)/(WTRAS(1,4)),'r'); %G13D - Mech 1
plot(GEFvals,100*WTRAS(:,5)/(WTRAS(1,5)),'m'); %G13D - Mech 2
set(gca, 'XDir','reverse');
set(gca,'XScale','log');
xlabel('(EGFR) SOS1/2 Activity');
ylabel('WT RAS-GTP (% max)');

figure(753)
hold on
plot(GEFvals,100*MUTRAS(:,1)/(MUTRAS(1,1)),'k');
plot(GEFvals,100*MUTRAS(:,2)/(MUTRAS(1,2)),'g');
plot(GEFvals,100*MUTRAS(:,3)/(MUTRAS(1,3)),'b');
plot(GEFvals,100*MUTRAS(:,4)/(MUTRAS(1,4)),'r'); %G13D - Mech 1
plot(GEFvals,100*MUTRAS(:,5)/(MUTRAS(1,5)),'m'); %G13D - Mech 2
set(gca, 'XDir','reverse');
set(gca,'XScale','log');
xlabel('(EGFR) SOS1/2 Activity');
ylabel('MUTANT RAS-GTP (% max)');