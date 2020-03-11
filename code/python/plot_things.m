cla
hold on;
ks = 3:6;

for k = ks
A = csvread(['results/tmp_N100_k', num2str(k), '_L1000.txt']);
plot(A(:,1),A(:,2).^2,'.');
end
xlabel('edit distance');
ylabel('squared inversion distance');
legends(strsplit(num2str(ks)));
%%
A = csvread(['results/tmp_N100_k3_L120_D70.txt']);
figure;
plot(A(:,1),A(:,3)/max(A(:,3)),'.');
hold on;
plot(A(:,1),A(:,4).^.5/max(A(:,4).^.5),'.');


%%
clf
t = 0:0.01:2;
for a =[2,4]
    y = 2/a/sqrt(2*pi).*t.^(1/a - 1).*exp(-1/2*t.^(2/a));
    plot(t,y);
    hold on;
end