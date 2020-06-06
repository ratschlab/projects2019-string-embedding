function sk = str_embed(x, siglen, Ind_tuples)
T = size(Ind_tuples,2);
dims = ones(1,T)*siglen;
sig_pow = siglen.^(0:(T-1));
sk = zeros(dims);
lim_inds = numel(x);
lim_inds = all(Ind_tuples<=lim_inds,2);
inds_pp = (x(Ind_tuples(lim_inds,:))-1)*sig_pow' + 1;
%     inds_pp = (x(Ind_tuples)-1)*sig_pow' + 1;
[uniq_inds,~,ic] = unique(inds_pp);
ind_cnts = accumarray(ic,1);
sk(uniq_inds) = ind_cnts;
