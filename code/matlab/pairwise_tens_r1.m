function dist_exact = pairwise_tens_r1(X,Alpha)
T = size(X{1},1);
D = size(X{1},2);
N = numel(X);

X_ext = zeros(N,D^T);
for i=1:N
    x = X{i};
    xext = x(1,:);
    for ti=2:T
        xext = xext(:)*x(ti,:);
    end
    X_ext(i,:) = xext(:);
end
dist_exact = zeros(N);
for i=1:N
    for j=i+1:N
        XD = X_ext(i,:)-X_ext(j,:);
        dist_exact(i,j) = sum(abs(XD).^Alpha).^(1/Alpha);
    end
end
dist_exact = dist_exact + dist_exact';
