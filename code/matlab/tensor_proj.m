function Pr = tensor_proj(D,T,M,Alpha,Disc,method)
dims = ones(1,T)*D;
dims_M = [dims, M];
pd = makedist('Stable','alpha',Alpha,'beta',0,'gam',1,'delta',0);
switch method
    case "naive"
        if Alpha==1
            Pr = randn(dims_M)./randn(dims_M);
        elseif Alpha==2
            Pr = randn(dims_M);
        end
    case "phideal"
        phase_ideal = rand(dims_M);
        Pr = icdf(pd,phase_ideal);
    case "phquasi"
        phase_quasi = zeros(dims);
        for j=1:M
            phase = rand_phase(D,T);
            if j==1
                phase_quasi = phase;
            else
                phase_quasi = cat(T+1,phase_quasi,phase);
            end
        end
        Pr = icdf(pd,phase_quasi);
    case "phdisc"
        phase_disc = zeros(dims);
        for j = 1:M
            res_phase = rand_phase(D,T,Disc);
            
            if j==1
                phase_disc = res_phase;
            else
                phase_disc = cat(T+1,phase_disc,res_phase);
            end
        end
        Pr = icdf(pd,phase_disc);
    otherwise
        error(['method ', method, ' does not exist']);
end
