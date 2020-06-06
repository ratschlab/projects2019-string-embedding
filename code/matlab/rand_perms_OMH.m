function Perms = rand_perms_OMH(M,L,siglen)

Perms = zeros(M,siglen,L);
for m=1:M
    Perms(m,:,:) = reshape(randperm(siglen*L),siglen,L);
end