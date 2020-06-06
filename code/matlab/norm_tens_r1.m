% it returns 1 scalar for T x D matrix, or N scalor for a cell 
function Norm = norm_tens_r1(X,Alpha)
if ~iscell(X) 
    X = {X};
end
N = numel(X);
Norm = zeros(1,N);
for i=1:N
    x = X{i};
    Norm(i) = prod(sum(abs(x).^Alpha,2)).^(1/Alpha); % compute norm 
end
