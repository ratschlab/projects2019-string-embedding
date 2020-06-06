function [opts] = opts_rp_sketch(vargs)

% dist_methods = {'median',"binary"};
is_dist_method = @(x) any(strcmpi(x,["median","binary","norm","hamming","LCS"]));
valid_alpha = @(x) isnumeric(x) && isscalar(x) && (x>0) ;

pars = inputParser;
addParameter(pars,'dist','norm',is_dist_method);
addParameter(pars,'alpha',2, valid_alpha);
parse(pars,vargs{:});
opts = pars.Results;
