% params 
clear;
N =200;
M = 200;
D = 100;
T = 2;
Disc = 8;
Alpha = 1;
Sim = .9985;   % adjacent similarity 

% Generate data
X = cell(1,N); % inputs
x_prev = (randn(T,D));
for i=1:N
    x = (randn(T,D)); % random  T x D
    x = (1-Sim) * x + Sim * x_prev; % smooth change  
    x = x/norm_tens_r1(x,Alpha)^(1/T); % normalize
    x_prev = x; % save previous 
    X{i} = x; % save sample i 
end

% initialize sketching
[Phase,iCDF] = rand_phase_disc(Disc,M,T,D,Alpha);
% Compute rank1 sketches
Sk = zeros(N,M);  
for i=1:N
    Sk(i,:) = sketch_rank1(X{i},Phase,iCDF);
end

% Exact & estimate rank1 dists
dist_exact = pairwise_tens_r1(X,Alpha);
dist_est = pairwise_dist(Sk,'dist','median');

%% visualize 
clf;
subplot(2,1,1);
plot_XY(dist_exact, dist_est,'Xlabel',"exact",'Ylabel',"Sk-est",'sort',1);
subplot(2,1,2);
plot_terrain(dist_exact, dist_est,'Xlabel',"exact",'Ylabel',"Sk-est");
