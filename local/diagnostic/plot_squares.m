function plot_squares(center_x, center_y, width, alpha, color)
% Plot squares on an existing figure.
%
%
%USAGE:
% PLOT_SQUARES(center_x, center_y, width, alpha)
%
%
%INPUT:
% center_x      /1xN vector/ x coordinate of marker center,
%                N is the number of squares
% center_y      /1xN vector/ y coordinate of marker center
% width         dimension of the marker /1xN vector/, if only a single 
%               integer is provided it will be used for every markers
%
%OPTIONAL:
% alpha         transparency of the marker's face color /[0,1] scalar/
% color         color of the marker's face fill /1x3 vector/ on [0,1]
%
%OUTPUT:
% Scatter plot with square markers and adjustable transparency.
% Transparency option is not available in _plot_ or _scatter_ plots.
% Due to its intended use only working as addition to existing figure.
% This plot preserves the axis limits of the original plot!
% Intended for plotting Monte Carlo simulation results.
%
% NOTES:
%  (1) The marker shape could be easily changed by adding verticies.
%  (2)_center_x_, _center_y_ and _width_ are reshaped to row vectors.
%

center_x = center_x(:).';
center_y = center_y(:).';
width    = width(:).';

if numel(center_x) ~= numel(center_y)
    error('''center_x'' and ''center_y'' should have the same number of elements!')
end

nsquare = size(center_x,2);
if isscalar(width)
    width = repmat(width,1,nsquare);
elseif numel(center_x) ~= numel(width)
    error('''width'' and ''center_x/y'' should have the same number of elements if ''width'' is not a scalar!')
end

if nargin < 5
    color = [0.7, 0.7, 0.7];
end
if nargin < 4
    alpha = 0.2;
end

fPos = get(gcf, 'Position');
% WARNING! scaling of the markers is based on the width and height of the 
% plot it is also dependent on _subplot_; the code is set to
% square arrangement

% need width, height in data values
xl = xlim();
yl = ylim();

w = width*(xl(2)-xl(1))/fPos(3);
h = width*(yl(2)-yl(1))/fPos(4);

% vertices
X = [center_x-w/2; center_x-w/2; center_x+w/2; center_x+w/2];
Y = [center_y-h/2; center_y+h/2; center_y+h/2; center_y-h/2];

% dummy
C = ones(size(X));

%plot patches
patch(X,Y,C,...
    'EdgeColor', 'none',...
    'FaceColor', color,...
    'FaceAlpha', alpha)

% set axis limits back to original!
% the plot will not adjust to points outside of the original
% plotting area!
xlim(xl)
ylim(yl)

end
