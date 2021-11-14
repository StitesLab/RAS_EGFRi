GEFvals=[10:-.25:1];
[m,n]=size(GEFvals);

RasIndices=[0,1,2,3,4,5,6,7,8,9];

for ii=1:n
    GEFfact=GEFvals(ii);

    for jj=1:10
    
        mutflag=RasIndices(jj);
        if jj==1;
            [KnownRasGTP(ii,jj),KnownEffRasGTP(ii,jj),KnownWTRasGTP(ii,jj),KnownMutRasGTP(ii,jj)]=ssRas(mutflag,0,1,GEFfact);
        else
            [KnownRasGTP(ii,jj),KnownEffRasGTP(ii,jj),KnownWTRasGTP(ii,jj),KnownMutRasGTP(ii,jj)]=ssRas(mutflag,.25,.75,GEFfact);
        end
    end
end


% at GEF=2x, index=33
KnownRasGTP=KnownRasGTP*100;
KnownRasGTPsub=KnownRasGTP(:,:)-(KnownRasGTP(37,:));

KnownEffRasGTP=KnownEffRasGTP*100;
KnownEffRasGTPsub=KnownEffRasGTP(:,:)-(KnownEffRasGTP(37,:));

KnownWTRasGTP=KnownWTRasGTP*100;
KnownWTRasGTPsub=KnownWTRasGTP(:,:)-(KnownWTRasGTP(37,:));

KnownMutRasGTP=KnownMutRasGTP*100;
KnownMutRasGTPsub=KnownMutRasGTP(:,:)-(KnownMutRasGTP(37,:));

KnownRasGTPfullnorm=KnownRasGTPsub(:,:)./(KnownRasGTPsub(1,:));
KnownRasGTPnorm=KnownRasGTP(:,:)./(KnownRasGTP(1,:));

KnownEffRasGTPfullnorm=KnownEffRasGTPsub(:,:)./(KnownEffRasGTPsub(1,:));
KnownEffRasGTPnorm=KnownEffRasGTP(:,:)./(KnownEffRasGTP(1,:));

KnownWTRasGTPfullnorm=KnownWTRasGTPsub(:,:)./(KnownWTRasGTPsub(1,:));
KnownWTRasGTPnorm=KnownWTRasGTP(:,:)./(KnownWTRasGTP(1,:));

KnownMutRasGTPfullnorm=KnownMutRasGTPsub(:,:)./(KnownMutRasGTPsub(1,:));
KnownMutRasGTPnorm=KnownMutRasGTP(:,:)./(KnownMutRasGTP(1,:));


Known2xRasGTP=KnownRasGTP(33,:);
Known2xEffRasGTP=KnownEffRasGTP(33,:);
Known2xWTRasGTP=KnownWTRasGTP(33,:);
Known2xMutRasGTP=KnownMutRasGTP(33,:);

Known2xRasGTPfullnorm=KnownRasGTPfullnorm(33,:);
Known2xRasGTPnorm=KnownRasGTPnorm(33,:);

Known2xEffRasGTPfullnorm=KnownEffRasGTPfullnorm(33,:);
Known2xEffRasGTPnorm=KnownEffRasGTPnorm(33,:);

Known2xWTRasGTPfullnorm=KnownWTRasGTPfullnorm(33,:);
Known2xWTRasGTPnorm=KnownWTRasGTPnorm(33,:);

Known2xMutRasGTPfullnorm=KnownMutRasGTPfullnorm(33,:);
Known2xMutRasGTPnorm=KnownMutRasGTPnorm(33,:);


Integral_KnownWTRasGTPfullnorm=sum(KnownWTRasGTPfullnorm);
Integral_KnownEffRasGTPfullnorm=sum(KnownEffRasGTPfullnorm);

Raslabels={'WT', 'G12V', 'G12D', 'G13D','Q61H','Q61K','Q61L','Q61P','Q61R','Q61W'};

G13Dintegral=Integral_KnownWTRasGTPfullnorm(4);
G13DintegralEff=Integral_KnownEffRasGTPfullnorm(4);

Diff(1)=Integral_KnownWTRasGTPfullnorm(2)-G13Dintegral;
Diff(2)=Integral_KnownWTRasGTPfullnorm(3)-G13Dintegral;

DiffEff(1)=Integral_KnownEffRasGTPfullnorm(2)-G13DintegralEff;
DiffEff(2)=Integral_KnownEffRasGTPfullnorm(3)-G13DintegralEff;

Diff2x(1)=Known2xWTRasGTPnorm(2)-Known2xWTRasGTPnorm(4);
Diff2x(2)=Known2xWTRasGTPnorm(3)-Known2xWTRasGTPnorm(4);
Diff2xfullnorm(1)=Known2xWTRasGTPfullnorm(2)-Known2xWTRasGTPfullnorm(4);
Diff2xfullnorm(2)=Known2xWTRasGTPfullnorm(3)-Known2xWTRasGTPfullnorm(4);


%%% Threshold for sensitiviyt - allows Q61K and Q61R but not others

MaxSens=min(Diff)/10;

MaxSensEff=min(DiffEff)/10;

MaxDiff2x=min(Diff2x)/10;
MaxDiff2xfullnorm=min(Diff2xfullnorm)/10;

%%%
%%% Find those that are sensitive

load DRs_RAND_MUTS_output_1e2

RasGTPsub=RasGTP(:,:)-(RasGTP(37,:));
RasEffsub=RasEff(:,:)-(RasEff(37,:));
WTRasGTPsub=WTRasGTP(:,:)-WTRasGTP(37,:);
MutRasGTPsub=MutRasGTP(:,:)-MutRasGTP(37,:);

RasGTPfullnorm=RasGTPsub(:,:)./(RasGTPsub(1,:));
RasEfffullnorm=RasEffsub(:,:)./(RasEffsub(1,:));
WTRasGTPfullnorm=WTRasGTPsub(:,:)./WTRasGTPsub(1,:);
MutRasGTPfullnorm=MutRasGTPsub(:,:)./MutRasGTPsub(1,:);

RasGTPnorm=RasGTP(:,:)./(RasGTP(1,:));
RasEffnorm=RasEff(:,:)./(RasEff(1,:));
WTRasGTPnorm=WTRasGTP(:,:)./WTRasGTP(1,:);
MutRasGTPnorm=MutRasGTP(:,:)./MutRasGTP(1,:);


Integral_WTRasGTPfullnorm=sum(WTRasGTPfullnorm);

Diff_Integral_WTRas=Integral_WTRasGTPfullnorm-G13Dintegral;

%%%%%%
zz=find(Diff_Integral_WTRas<MaxSens);

Sens_Rand_A=RANGE_RAND(zz,:);
Sens_WTDRs_A=WTRasGTPnorm(:,zz);

Sens_RasGTP_A=RasGTP(:,zz);
Sens_RasEff_A=RasEff(:,zz);
Sens_WTRasGTP_A=WTRasGTP(:,zz);
Sens_MutRasGTP_A=MutRasGTP(:,zz);

Diff_WT2x(:)=Sens_WTDRs_A(33,:)-Known2xWTRasGTPnorm(4);
yy=find(Diff_WT2x<MaxDiff2x);

Sens_Rand=Sens_Rand_A(yy,:);

Sens_RasGTP=Sens_RasGTP_A(:,yy);
Sens_RasEff=Sens_RasEff_A(:,yy);
Sens_WTRasGTP=Sens_WTRasGTP_A(:,yy);
Sens_MutRasGTP=Sens_MutRasGTP_A(:,yy);

Kms=Sens_Rand(:,7);
smallKms=find(Kms<1);
Kmvals=Sens_Rand(smallKms,7);

Sens_Rand2=Sens_Rand(smallKms,:);

%%% Calculate Kd values
MutKds=Sens_Rand2(:,9)./Sens_Rand2(:,8);
MutKdsRel=MutKds/(80e-9);

save Sens_Rand_1e2.txt Sens_Rand -ASCII -double
save Sens_Rand2_1e2.txt Sens_Rand2 -ASCII -double

%figure(1003), semilogy(Sens_Rand','k.')
%figure(1004), semilogy(Sens_Rand2','k.')



