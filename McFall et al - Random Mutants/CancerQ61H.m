mkint=(1.83/20.93)*kint;
mkdissD=(9.6/6.5)*kdissD;
mkdissT=(9.2/15)*kdissT;
mkassD=kassD;
mkassT=kassT;
mkcat=mkint;
mKm=(25/35)*Km;
mkdissEff=kdissEff;
mkassEff=kassEff;
mH=(mkassD*mkdissT)/(mkdissD*mkassT);
mkD=3.9;
mKmD=3.86e-4/volscale;
mKmT=3e-4/volscale;
mkT=mkD*mKmT*mH/mKmD;