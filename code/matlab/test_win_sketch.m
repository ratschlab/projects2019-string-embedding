% params 
clear;
N =5;
L = 5; % length of string 
Alpha = 1;
mutate_rate = .6;   % adjacent similarity 
siglen = 10; % number of characters

% Generate data
X = cell(1,N); % inputs
x_prev = randi(siglen,1,L); % random  T x D
S_con = zeros(L,N);
for i=1:N
    x = x_prev; 
    mutate_inds = rand(1,L)<mutate_rate;
    num_mut = sum(mutate_inds);
    x(mutate_inds) = randi(siglen,1,num_mut);
    x_prev = x; % save previous 
    X{i} = x; % save sample i 
    S_con(:,i) = x;
end
S_con = S_con(:);


LCS_dist = pairwise_dist(X,'dist',"LCS");
hamming_dist = pairwise_dist(X,"dist","hamming");


%% Efficient sketching 
M = 1;
Disc = 5;
T = 3;

[Phase,iCDF] = rand_phase_disc(Disc,M,T,siglen,Alpha);
Sk = zeros(N,M);  
Counts = cell(1,N);
for i=1:N
    [Sk(i,:),Counts{i}] = sketch_str_tup(X{i},Phase,iCDF);
end
Sk_dist = pairwise_dist(Sk,'dist','median');

Sk2 = sketch_str_win(S_con,L,Phase,iCDF);
Sk2 = Sk2(L:L:end,:);
Sk_dist2 = pairwise_dist(Sk2,'dist','median');


clf
subplot 211;
plot_terrain(LCS_dist,Sk_dist);
subplot 212;
plot_terrain(LCS_dist,Sk_dist2, 'opacityGray',.015, 'sort',1 );
