function [Phase,iCDF,Bins] = rand_phase_disc(Disc,M,T,D,Alpha,nbins)
PDF = makedist('Stable','alpha',Alpha,'beta',0,'gam',1,'delta',0);
Phase = randi(Disc,M,T,D);
Bias = sum(rand(T,M),1);
Phase_uvals = mod((0:Disc-1)'+Bias,Disc)/Disc;
iCDF = icdf(PDF,Phase_uvals);
Bins = icdf(PDF,(0:nbins)/nbins);