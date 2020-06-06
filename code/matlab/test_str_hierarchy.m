% params
clear;
N =200;
L = 200; % length of string
Alpha = 1;
mutate_rate = .02;   % adjacent similarity
block_mute = 0.1;
max_blocks = 4;
siglen = 2; % number of characters


X = gen_strings(siglen,N,L,mutate_rate,block_mute,max_blocks);


LCS_dist = pairwise_dist(X,'dist',"LCS");

%% random projection sketching
M = 200;
nbins = 16;
nbins2 = 1;
Disc = 8;
T = 3;
Div = 5; % divide to `Div` sub-strings
OL = .5; % overlap with each succeeding window

W = round(L/Div);
Inc = round(W*(1-OL));
numInc = floor((L-W+1)/Inc);

[Phase,iCDF,Bins] = rand_phase_disc(Disc,M,T,siglen,Alpha,nbins);
[Phase2,iCDF2,Bins2] = rand_phase_disc(Disc,M,T,nbins,Alpha,nbins2);

Sk = cell(N,1);
for xi=1:N
    x = X{xi};
    Sk1 = zeros(M,numInc);
    for i=1:numInc
        x_sub = x((1:W)+(i-1)*Inc);
        Sk1(:,i) = sketch_str_tup(x_sub,Phase,iCDF,Bins);
    end
    Sk2 = zeros(M,1);
    for mi=1:M
        subSk = Sk1(mi,:);
        Sk2(mi,:) = sketch_str_tup(Sk1(mi,:),Phase2(mi,:,:),iCDF2(:,mi),[]);
    end
    Sk{xi} = Sk2(:);
end
Sk_dist = pairwise_dist(Sk,'dist','binary');

subplot 221;
cla;
plot_XY(LCS_dist,Sk_dist);
subplot 222;
cla;
plot_terrain(LCS_dist,Sk_dist, 'opacityGray',.015, 'sort',1 );

%% Ordered MinHash
Tupple = 1;

Perms = rand_perms_OMH(M,L,siglen);
Sk_OMH = cell(N,1);
for i=1:N
    Sk_OMH{i} = sketch_str_OMH(X{i},Perms,T,Tupple);
end
OMH_dist = pairwise_dist(Sk_OMH,'dist',"hamming");


subplot 223;
cla;
plot_XY(LCS_dist,OMH_dist);
subplot 224;
cla;
plot_terrain(LCS_dist,OMH_dist, 'opacityGray',.015, 'sort',1 );

disp('rand proj & OMH, spearman')
corr(LCS_dist(:),Sk_dist(:),'type','Spearman')
corr(LCS_dist(:),OMH_dist(:),'type','Spearman') 
disp('rand proj & OMH, Pearson')
corr(LCS_dist(:),Sk_dist(:),'type','Pearson')
corr(LCS_dist(:),OMH_dist(:),'type','Pearson')

% disp(['time for random projection = ', num2str(rptime)]);
% disp(['time for ordered min-hash = ', num2str(omhtime)]);
