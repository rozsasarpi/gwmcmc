function q = quan_tile(X, p, method)
% Quantiles of a data set
%
%USAGE:
% q = QUAN_TILE(x, p, method)
%
%INPUT:
% x         sample vector /vector/
% p         probabilities of required quantiles [0,1] /vector, 1xN/
%
%OPTIONAL:
% method    Plotting position, which formula is used for constructing the 
%           empirical cumulative distribution function
%   1       empiF = (i-2/5)/(n+1/5) [default]
%   2       empiF = (i-0.5)/(n); same as in Matlab's Statistical Toolbox _quantile.m_
%
%OUTPUT:
% q         quantiles /vector, 1xN/
%
%NOTES:
%   Quantiles are specified using cumulative probabilities, from 0 to 1.
%   For an N element vector _x_, QUANTILE computes quantiles as follows:
%      1) Linear interpolation is used to compute quantiles for
%         probabilities between _min(empiF)_ and _max(empiF)_
%      2) The minimum or maximum values in _x_ are assigned to quantiles
%         for probabilities outside that range.
%      3) NaN and Inf values are omitted in the calculations.
%
%   Simple implementation to remove the dependence on toolboxes.

if nargin < 3
    method = 1;
end

if any(p > 1) || any(p < 0)
    error('''p'' values should be in interval [0,1]!')
end

if ~ismatrix(X) || isempty(X)
    error('''x'' should be a vector or a matrix!')
end
if ~isvector(p) || isempty(p)
    error('''p'' should be a vector!')
end

% discard NaN and Inf values from calculation
idx = isfinite(X);
if any(~idx)
    X(~idx) = [];
end

n   = size(X,1);
sX  = sort(X);

switch method
    case 1
        empiF = ((1:n)-2/5)/(n+1/5);        
    case 2
        empiF = (0.5:(n-0.5))/(n);
    otherwise
       error(['Unknown method identifier: ', num2str(method)]) 
end

idx_below       = p < min(empiF);
idx_above       = p > max(empiF);
idx_in          = not(idx_below | idx_above);

q               = nan(length(p),size(X,2));
% this ugly _if_ is needed to protect against dimension mismatch
% in case of empty-matrix..
if any(idx_below)
    q(idx_below,:)  = min(X);
end
if any(idx_above)
    q(idx_above,:)  = max(X);
end
if any(idx_in)
    q(idx_in,:)     = interp1(empiF, sX, p(idx_in));
end

end
