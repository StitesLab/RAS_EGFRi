load RANGE_RAND


GEFvals=[10:-.25:1];
[m,n]=size(GEFvals);

[aaa,bbb]=size(RANGE_RAND);
No_Muts=max(aaa,bbb);

RasGTP=zeros(n,No_Muts);
RasEff=RasGTP;
WTRasGTP=RasGTP;
MutRasGTP=RasGTP;


for ii=1:n;
    parfor jj=1:No_Muts
        [RasGTP(ii,jj),RasEff(ii,jj),WTRasGTP(ii,jj),MutRasGTP(ii,jj)]=ssRas_RAND(RANGE_RAND(jj,:),0.25,0.75,GEFvals(ii));
    end
end

RasGTP=RasGTP*100;
RasEff=RasEff*100;
WTRasGTP=WTRasGTP*100;
MutRasGTP=MutRasGTP*100;

save DRs_RAND_MUTS_output
