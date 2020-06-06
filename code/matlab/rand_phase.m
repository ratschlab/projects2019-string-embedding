function res_phase = rand_phase(d,t,disc)

is_real = false;
if nargin < 3
    is_real = true;
    disc = 1;
end

inds = ones(1,t);
inds(1) = d;
if is_real
    phase = rand(inds) ;
else
    phase = randi(disc,inds) ;
end

for k=2:t
    inds = ones(1,t);
    inds(k) = d;
    if is_real
        phase1D = rand(inds) ;
    else
        phase1D = randi(disc,inds) ;
    end
    phase = phase + phase1D;
end
phase_bias = 0; 
if ~is_real
    phase_bias = sum(rand(1,t));
end
res_phase = mod(mod(phase-t,disc)+phase_bias,disc)/disc;