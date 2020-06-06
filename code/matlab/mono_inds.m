function Ind_tuples = mono_inds(L,T)

num_tup = nchoosek(L,T);
Ind_tuples = zeros(num_tup,T);
tup = ones(1,T+1);
i = 1;
while (tup(T+1)==1) 
    % check if monotonic
    if all(diff(tup(1:T))>0)
        Ind_tuples(i,:) = tup(1:T);
        i = i + 1;
    end
    % plus one index 
    tup(1) = tup(1) + 1;
    ii = 1;
    while tup(ii)>L
        tup(ii) = 1;
        ii = ii + 1;
        tup(ii) = tup(ii) + 1;
    end
end
