function rho = correlation(x,y,type)
% Linear or rank correlation between two vectors
%
%USAGE:
%
% rho = CORRELATION(x,y,type)
%
%INPUT:
% x
% y
%
%OPTIONAL:
% type      [default='spearman']
%
%NOTES:
% Simple implementation to remove the dependence on toolboxes...
%

if nargin < 3
    type = 'spearman';
end

% reshape _x_ and _y_ into column vectors
x = x(:);
y = y(:);

nx = numel(x);
ny = numel(y);
if nx ~= ny
    error('x and y should have the same size!')
end

switch lower(type)
    case {'spearman', 'spear', 'spearm', 's'}
        % does not handle ties!
        [~, ix] = sort(x);
        [~, rx] = sort(ix);
        [~, iy] = sort(y);
        [~, ry] = sort(iy);
        d       = rx - ry;
        n       = numel(x);
        
        rho     = 1-6*sum(d.^2)/(n*(n^2-1));
%     case {'kendall', 'ken', 'k'}
%         
    case {'pearson', 'pears', 'p'}
        cent_x = x - mean(x);
        cent_y = y - mean(y);
        
        rho = cent_x.'*cent_y/(norm(cent_x)*norm(cent_y));
       
    otherwise
        error(['Unknown type: ', type])
end