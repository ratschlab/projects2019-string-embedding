function D = pairwise_dist(X,varargin)
opts = opts_rp_sketch(varargin);
a = opts.alpha;
N = size(X,1);
D = zeros(N,N);
for i=1:N
    xi = getXi(X,i);
    for j=i+1:N
        xj = getXi(X,j);
        switch opts.dist
            case 'median'
                D(i,j)= median(abs(xi-xj));
            case 'binary'
                D(i,j)= sum(abs((xi>0)-(xj>0)));
            case 'norm'
                D(i,j) = sum(abs(xi-xj).^a,'all').^(1/a);
            case 'hamming'
                D(i,j) = sum(xi~=xj,'all');
            case 'LCS'
                D(i,j) = 1- LCS(xi,xj);
        end
    end
end

D = D + D';