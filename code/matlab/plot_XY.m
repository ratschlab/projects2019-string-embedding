function  plot_XY(varargin)
[X,Y,opts] = opts_plotXY(varargin);

n = size(X,1);
if opts.sort
    [~,I] = sort(X,2);
    for i=1:n
        Y(i,:) = Y(i,I(i,:));
    end
end
if opts.logy
    Y = log(Y);
end
imagesc(Y);
xlabel(opts.Xlabel);
ylabel(opts.Ylabel);
if opts.cbar
    colorbar;
end
colormap(opts.cmap);
