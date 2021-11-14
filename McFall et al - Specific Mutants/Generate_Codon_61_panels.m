GEFvals=[10:-.25:1];
[m,n]=size(GEFvals);

RasIndices=[0,1,2,3,4,5,6,7,8,9,101,102];

for ii=1:n
    GEFfact=GEFvals(ii);

    for jj=1:12
    
        mutflag=RasIndices(jj);
        if jj==1;
            [RasGTP(ii,jj),RasEff(ii,jj),WTRasGTP(ii,jj),MutRasGTP(ii,jj)]=ssRas(mutflag,0,1,GEFfact);
        else
            [RasGTP(ii,jj),RasEff(ii,jj),WTRasGTP(ii,jj),MutRasGTP(ii,jj)]=ssRas(mutflag,.25,.75,GEFfact);
        end
    end
end

RasGTP=RasGTP*100;
RasEff=RasEff*100;
WTRasGTP=WTRasGTP*100;
MutRasGTP=MutRasGTP*100;

Raslabels={'WT', 'G12V', 'G12D', 'G13D','Q61H','Q61K','Q61L','Q61P','Q61R','Q61W','Q61K/WT','Q61R/WT'};

stimlabels={'low','high'};

figure(10);
for ii=5:10
    subplot(3,2,ii-4);
    hold on
plot(GEFvals,100*RasGTP(:,2)/RasGTP(1,2),'k'); % G12V
plot(GEFvals,100*RasGTP(:,4)/RasGTP(1,4),'b'); %G13D
plot(GEFvals,100*RasGTP(:,ii)/RasGTP(1,ii),'r');

axis([1 10 25 100]);
set(gca,'xtick',[1,10],'xticklabel',stimlabels)
set(gca, 'XDir','reverse','Xscale','log')
xlabel('(EGFR) SOS1/2 activity');
ylabel('RasGTP (% total Ras)');
end


figure(11);
for ii=5:10
    subplot(3,2,ii-4);
    hold on
plot(GEFvals,100*MutRasGTP(:,2)/MutRasGTP(1,2),'k');
plot(GEFvals,100*MutRasGTP(:,4)/MutRasGTP(1,4),'b'); %G13D
plot(GEFvals,100*MutRasGTP(:,ii)/MutRasGTP(1,ii),'r');

set(gca,'xtick',[1,10],'xticklabel',stimlabels)
set(gca, 'XDir','reverse','Xscale','log')
axis([1 10 0 120]);
xlabel('(EGFR) SOS1/2 activity');
ylabel('Mutant RasGTP (% max)');
end

figure(12);
for ii=5:10;
    subplot(3,2,ii-4);
    hold on
plot(GEFvals,100*WTRasGTP(:,2)/WTRasGTP(1,2),'k');
plot(GEFvals,100*WTRasGTP(:,4)/WTRasGTP(1,4),'b'); %G13D
plot(GEFvals,100*WTRasGTP(:,ii)/WTRasGTP(1,ii),'r');

set(gca,'xtick',[1,10],'xticklabel',stimlabels)
set(gca, 'XDir','reverse','Xscale','log')
axis([1 10 0 100]);
xlabel('(EGFR) SOS1/2 activity');
ylabel('WT RasGTP (% max)');
end

figure(61); 
subplot(2,1,1)
% Dose Response - Q61K with WT NF1 affinity
hold on
plot(GEFvals,100*RasGTP(:,2)/RasGTP(1,2),'k'); %G12V
plot(GEFvals,100*RasGTP(:,4)/RasGTP(1,4),'b'); %G13D
plot(GEFvals,100*RasGTP(:,6)/RasGTP(1,6),'r'); %Q61K
plot(GEFvals,100*RasGTP(:,11)/RasGTP(1,11),'r:'); %Q61K with WT NF1 affinity
axis([1 10 0 100]);
set(gca,'xtick',[1,10],'xticklabel',stimlabels)
set(gca, 'XDir','reverse','Xscale','log')
xlabel('(EGFR) SOS1 activity');
ylabel('RasGTP (% max)');

subplot(2,1,2)
hold on
plot(GEFvals,100*RasGTP(:,2)/RasGTP(1,2),'k'); %G12V
plot(GEFvals,100*RasGTP(:,4)/RasGTP(1,4),'b'); %G13D
plot(GEFvals,100*RasGTP(:,9)/RasGTP(1,9),'r'); %Q61R
plot(GEFvals,100*RasGTP(:,12)/RasGTP(1,12),'r:'); %Q61R with WT NF1 affinity
axis([1 10 0 100]);
set(gca,'xtick',[1,10],'xticklabel',stimlabels)
set(gca, 'XDir','reverse','Xscale','log')
xlabel('(EGFR) SOS1 activity');
ylabel('RasGTP (% max)');