function  plot_terrain(varargin)
[X,Y,opts] = opts_plotXY(varargin);
assert(size(Y,1)==size(Y,1), 'X and Y must match');

n = size(X,1);
[sX,I] = sort(X);
sY = Y;
for i=1:n
    sY(:,i) = Y(I(:,i),i);
end

hold on;
rcols = randi(n,1,opts.numcol);
rcols = sort(rcols);
for i=1:opts.numcol
    ri = rcols(i);
    plot(sX(:,ri),sY(:,ri),'.','Color',[rand(1,3),opts.opacityColor],'markersize',8,'LineWidth',1.5,'LineStyle','-');
end
plot(sX,sY,'Color',[0,0,0,opts.opacityGray]);
xlabel(opts.Xlabel);
ylabel(opts.Ylabel);

legend(strsplit(num2str(rcols)),'Location','northwest');
hold off;

