clear;
N =100;
M = 50;
D = 25;
T = 3;
Disc = 5;
Alpha = 1;
adj_sim = .95;   % adjacent similarity
scale = rand(1,N);

proj_methods = ["Naive","phIdeal", "phQuasi","phDisc"];
num_methods = numel(proj_methods);
Projs = cell(1,num_methods);

for method_j=1:num_methods
    method = lower(proj_methods(method_j));
    Pr = tensor_proj(D,T,M,Alpha,Disc,method);
    Projs{method_j} = Pr;
    
end

dims = ones(1,T)*D;
X = cell(N,1);
for i=1:N
    x = scale(i)*rand(dims);
    x = x/sum(abs(x).^Alpha,'all').^(1/Alpha);
    if i==1
        xp = x;
    else
        x = adj_sim * xp + (1-adj_sim) * x;
        xp = x;
    end
    X{i} = x;
end

copy_num = ones(1,T+1);
copy_num(T+1) = M;
Sk = zeros(num_methods,N,M);
for i=1:N
    x = X{i};
    xrep = repmat(x,copy_num);
    for method_j=1:num_methods
        Sk(method_j,i,:) = sum(xrep.*Projs{method_j},1:T);
    end
end


dist_est = cell(1,num_methods);
for method_j=1:num_methods
    dist_est{method_j} = pairwise_dist(squeeze(Sk(method_j,:,:)),'dist','median');
end

dist_exact = pairwise_dist(X,'dist',"norm",'alpha',Alpha);

%%
clf
hold on;
for method_j = 1:num_methods
    subplot(2,2,method_j);
    
    name = proj_methods(method_j);
    D2 = dist_est{method_j};
    plot_terrain(dist_exact,D2, "xlabel","exact","ylabel",name);
    Corr = corr(dist_exact(:),D2(:));
    title(['Corr = ', num2str(Corr)]);
end

