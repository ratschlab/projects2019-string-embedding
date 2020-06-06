function sk = sketch_rank1(x,Phase,iCDF)
Disc = size(iCDF,1);
M = size(iCDF,2);
T = size(Phase,2);

Coefs = zeros(Disc,M); % Disc x M
for mi=1:M
    phase_mi = squeeze(Phase(mi,:,:)); % phase2: t x d 
    coefs_mi = zeros(Disc,T);
    for di=1:Disc
        coefs_mi(di,:) = sum(x.*(phase_mi==di),2);
    end
    fft_coef = fft(coefs_mi);
    fft_prod = prod(fft_coef,2);
    Pr = ifft(fft_prod);
    Coefs(:,mi) = Pr;
end
sk = sum(Coefs.*iCDF,1);