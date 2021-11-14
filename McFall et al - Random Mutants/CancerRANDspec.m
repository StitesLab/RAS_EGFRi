    mkint=RANDfactor(1)*kint;

    mkdissD=RANDfactor(2)*kdissD;

    mkdissT=RANDfactor(3)*kdissT;

    mkassD=RANDfactor(4)*kassD;

    mkassT=RANDfactor(5)*kassT;

    mkcat=RANDfactor(6)*kcat;

    mKm=RANDfactor(7)*Km;

    mkassEff=RANDfactor(8)*kassEff;

    mkdissEff=RANDfactor(9)*kdissEff;

    mkD=RANDfactor(10)*kD;

    mKmD=RANDfactor(11)*KmD;

    mKmT=RANDfactor(12)*KmT;

mHaldaneint=(mkassD*mkdissT)/(mkdissD*mkassT);
mkT=mkD*mKmT*mHaldaneint/mKmD;