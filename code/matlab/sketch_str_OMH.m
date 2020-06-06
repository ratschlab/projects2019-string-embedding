function sk = sketch_str_OMH(x,Perms,T,Tup)
if nargin<4
    Tup = 1;
end

% Perms = M x siglen x L
[M,siglen,L] = size(Perms);

l = length(x);
xc = zeros(1,l);
count = ones(1,siglen);
sig2ind = zeros(siglen,L);
for i=1:l
    a = x(i);
    xc(i) = count(a);
    sig2ind(a,xc(i)) = i;
    count(a) = count(a) + 1;
end
sk = zeros(T,M);
for mi=1:M
    inds = sub2ind(size(Perms),mi*ones(1,l),x,xc);
    P = Perms(inds);
    [~,sI] = sort(P);
    inds = inds(sI(1:T));
    [~,xL,xcL] = ind2sub(size(Perms),inds);
    indL = sig2ind(sub2ind(size(sig2ind),xL,xcL));
    indL = sort(indL);
    sk(:,mi) = x(indL);
end
if Tup==1
    sk = sk(:);
else
    Mask = siglen.^(0:T-1);
    sk = Mask*sk;
end

