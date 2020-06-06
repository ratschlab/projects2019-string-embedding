function Sk = sketch_str_win(x,W, Phase,iCDF)
L = length(x);
Disc = size(iCDF,1);
M = size(iCDF,2);
T = size(Phase,2);

Count = zeros(Disc,M,T); % Disc x M x T
DM = Disc*M;
[Dind,Mind] = ind2sub([Disc,M],1:DM);
Tind = ones(1,DM);
sqPhase = squeeze(Phase);
Sk = zeros(L,M);
XP = zeros(T,length(x));
for ti=1:T
    XP(ti,:) = sqPhase(sub2ind(size(sqPhase),ti*ones(1,length(x)),x'));
end
Counts = cell(L,1);
for i=1:L
    % Count[i][ph + ph(x_i) \mod Disc ]+= Cout[i-1][ ph ] , Disc)
    for t=T:-1:2
        xi_ph = squeeze(Phase(:,t,x(i)));
        xi_ph = repmat(xi_ph',Disc,1);
        Dind_sh = mod((Dind(:)-1) + (xi_ph(:)-1),Disc)+1;
        iiH = sub2ind(size(Count),Dind_sh(:),Mind(:),t*Tind(:));
        iiL = sub2ind(size(Count),Dind(:),Mind(:),(t-1)*Tind(:));
        Count(iiH) = Count(iiH) + Count(iiL);
    end
    iiL = sub2ind(size(Count),squeeze(Phase(:,1,x(i)))', 1:M, 1*ones(1,M));
    Count(iiL) = Count(iiL) + 1;
    
    
    % Count[i] -= mod( Count
    if i>W
        j = i - W;
        jjL = sub2ind(size(Count),squeeze(Phase(:,1,x(j)))', 1:M, 1*ones(1,M));
        Count(jjL) = Count(jjL) - 1;
        for t=2:T
            xj_ph = squeeze(Phase(:,t,x(j)));
            xj_ph = repmat(xj_ph',Disc,1);
            Dind_sh = mod((Dind(:)-1) + (xj_ph(:)-1),Disc)+1;
            jjL = sub2ind(size(Count),Dind(:),Mind(:),(t-1)*Tind(:));
            jjH = sub2ind(size(Count),Dind_sh(:),Mind(:),t*Tind(:));
            Count(jjH) = Count(jjH) - Count(jjL);
            
        end
    end
    % save the sketch
    CountT = squeeze(Count(:,:,T));
    Counts{i} = CountT;
    Sk(i,:) = sum(CountT.*iCDF);
end