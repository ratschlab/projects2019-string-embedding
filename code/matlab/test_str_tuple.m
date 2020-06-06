% params
clear;
N =200;
L = 200; % length of string
Alpha = 1;
mutate_rate = .02;   % adjacent similarity
block_mute = 0.05;
max_blocks = 1;
siglen = 4; % number of characters


X = gen_strings(siglen,N,L,mutate_rate,block_mute,max_blocks);


LCS_dist = pairwise_dist(X,'dist',"LCS");

%% Embed in a tuple
T = 2;
Ind_tuples = mono_inds(L,T);
tup_embed = cell(N,1);
for i=1:N
    tup_embed{i} = str_embed(X{i}, siglen, Ind_tuples);
end
tup_dist = pairwise_dist(tup_embed,'dist',"norm","alpha",1);

plot_XY(LCS_dist,tup_dist);
plot_terrain(LCS_dist,tup_dist, 'opacityGray',.03, 'sort',1 );

%% random projection sketching
M = 100;
Disc = 8;
T = 4;
nbins = 1;

tic
[Phase,iCDF,Bins] = rand_phase_disc(Disc,M,T,siglen,Alpha,nbins);
Sk = zeros(N,M);
for i=1:N
    Sk(i,:) = sketch_str_tup(X{i},Phase,iCDF,Bins);
end
if nbins==1
    Sk_dist = pairwise_dist(Sk,'dist','median');
else
    Sk_dist = pairwise_dist(Sk,'dist','hamming');
end
rptime = toc;


clf
subplot 221;
plot_XY(LCS_dist,Sk_dist);
subplot 222;
plot_terrain(LCS_dist,Sk_dist, 'opacityGray',.015, 'sort',1 );

% Ordered MinHash
Tupple = 0;

tic
Perms = rand_perms_OMH(M,L,siglen);
Sk_OMH = cell(N,1);
for i=1:N
    Sk_OMH{i} = sketch_str_OMH(X{i},Perms,T,Tupple);
end
OMH_dist = pairwise_dist(Sk_OMH,'dist',"hamming");
omhtime = toc;


clf
subplot 221;
plot_XY(LCS_dist,Sk_dist);
subplot 222;
plot_terrain(LCS_dist,Sk_dist, 'opacityGray',.015, 'sort',1 );

subplot 223;
plot_XY(LCS_dist,OMH_dist);
subplot 224;
plot_terrain(LCS_dist,OMH_dist, 'opacityGray',.015, 'sort',1 );

disp('rand proj & OMH, spearman')
corr(LCS_dist(:),Sk_dist(:),'type','Spearman')
corr(LCS_dist(:),OMH_dist(:),'type','Spearman') 
disp('rand proj & OMH, Pearson')
corr(LCS_dist(:),Sk_dist(:),'type','Pearson')
corr(LCS_dist(:),OMH_dist(:),'type','Pearson')

disp(['time for random projection = ', num2str(rptime)]);
disp(['time for ordered min-hash = ', num2str(omhtime)]);
