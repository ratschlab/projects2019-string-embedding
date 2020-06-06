function sk = sketch_str_tup(x,Phase,iCDF, Bins)
L = length(x);
Disc = size(iCDF,1);
M = size(iCDF,2);
T = size(Phase,2);

Count = zeros(Disc,M,T); % Disc x M x T
DM = Disc*M;
[Dind,Mind] = ind2sub([Disc,M],1:DM);
Tind = ones(1,DM);
for i=1:L
    for t=T:-1:2
        xi_ph = squeeze(Phase(:,t,x(i)));
        xi_ph = repmat(xi_ph',Disc,1);
        iiH = sub2ind(size(Count),Dind(:),Mind(:),t*Tind(:));
        Dind_sh = mod((Dind(:)-1) + (xi_ph(:)-1),Disc)+1;
        iiL = sub2ind(size(Count),Dind_sh(:),Mind(:),(t-1)*Tind(:));
        Count(iiH) = Count(iiH) + Count(iiL);
    end
    iiL = sub2ind(size(Count),squeeze(Phase(:,1,x(i)))', 1:M, 1*ones(1,M));
    Count(iiL) = Count(iiL) + 1;
end
CountT = squeeze(Count(:,:,T));
sk = sum(CountT.*iCDF,1); 
sk = sk./sum(CountT,1); % normalize
if numel(Bins)>2
    sk = sum(sk>Bins');
end
