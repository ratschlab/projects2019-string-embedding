clear;
n =200;
m = 100;
d = 100;
t = 2;
disc = 5;


% gauss_proj = randn(d,d,m);
cauchy_proj = randn(d,d,m)./randn(d,d,m);
cauchy_phase = rand(d,d,m);
cauchy_quant = tan(pi*(cauchy_phase-.5));
cauchy_quasi_phase = zeros(d,d,m);
for j=1:m
    cauchy_quasi_phase(:,:,j) = mod(rand(1,d)+rand(d,1),1);
end
cauchy_quasi_quant = tan(pi*(cauchy_quasi_phase-.5));

cauchy_quasi_phase_disc = zeros(d,d,m);
for j = 1:m
    phase_disc = randi(disc,d,1)+randi(disc,1,d);
    phase_bias = rand+rand;
    cauchy_quasi_phase_disc(:,:,j) = mod(phase_disc+phase_bias,disc)/disc;  
end
cauchy_quasi_quant_disc = tan(pi*(cauchy_quasi_phase_disc-.5));

C = eye(d);
for i=1:d
    for j=1:d
        C(i,j) = exp(-abs(i-j));
    end
end

res = zeros(5,n);
X = zeros(n,d,d);
Y = zeros(5,n,m);
for i=1:n
    scale = rand;
    x = scale*randn(d,d)*C;
    X(i,:,:) = x;
    xrep = (repmat(x,[1,1,m]));
    Y(1,i,:) = sum(xrep.*cauchy_proj,1:2);
    Y(2,i,:) = sum(xrep.*cauchy_quant,1:2);
    Y(3,i,:) = sum(xrep.*cauchy_quasi_quant,1:2);
    Y(4,i,:) = sum(xrep.*cauchy_quasi_quant_disc,1:2);
    
    res(1,i) = sum(abs(x),'all');
    res(2,i) = median(abs(Y(1,i,:)));
    res(3,i) = median(abs(Y(2,i,:)));
    res(4,i) = median(abs(Y(3,i,:)));
    res(5,i) = median(abs(Y(4,i,:)));
end

%%
clf
hold on;
method_name = {'','naive','phase-ideal','phase-quasi',['phase-disc: ',num2str(disc)]};
for j = 2:size(res,1)
    subplot(3,2,j-1);
    plot(res(1,:),res(j,:),'.','MarkerSize',8);
    title(method_name{j});
end
for j = 2:size(res,1)
    corr(res(1,:)',res(j,:)')
end

