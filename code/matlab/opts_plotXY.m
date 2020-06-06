function [X,Y,opts] = opts_plotXY(vargs)

validBinary = @(x) isnumeric(x) && isscalar(x);
validSquare = @(x) isnumeric(x) && (size(x,1)==size(x,2));
validInteger = @(x) isinteger(x) && isscalar(x);
number01 = @(x) isnumeric(x) && isscalar(x) && (x>=0 && x<=1);

pars = inputParser;
addRequired(pars, 'X',validSquare);
addRequired(pars, 'Y',validSquare);
addParameter(pars,'Xlabel', 'X',@isstring);
addParameter(pars,'Ylabel', 'Y',@isstring);
addParameter(pars,'cbar', 0, @isnumeric);
addParameter(pars,'sort',1,validBinary);
addParameter(pars,'logx',0,validBinary);
addParameter(pars,'logy',0,validBinary);
addParameter(pars,'opacityColor',.3,number01);
addParameter(pars,'opacityGray',.01,number01);
addParameter(pars,'cmap','jet',@isstring);
addParameter(pars,'numcol', 5,validInteger);
parse(pars,vargs{:});
opts = pars.Results;
X = opts.X;
Y = opts.Y;
assert(size(Y,1)==size(Y,1), 'X and Y must match');
