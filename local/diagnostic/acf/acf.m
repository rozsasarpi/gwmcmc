function varargout = acf(y,nLags)
% ACF - Compute Autocorrelations Through nLags lags
%
%USAGE:
% [myacf, lags, bounds] = ACF(y,Lags)
%
%INPUT:
% y         series to compute acf for, /nx1 column vector
% nLags     total number of lags, 1x1 integer
%
%OUTPUT:
% myacf - px1 vector containing autocorrelations
%        (Lag 0 also computed)
%
%
% If no output argument is specified it creates a plot, similar to _autocorr_
% of Econometrics Toolbox
%
%EXAMPLE:
% acf(randn(100,1), 10)
%
% NOTE:
% Majority of the code is written by Calvin Price
% mainly the plot is modified to get similar outcode as _autocorr_
% http://www.mathworks.com/matlabcentral/fileexchange/30540-autocorrelation-function--acf-

%
%TODO: check the user claim from FE (_autocorr_ gives the same results, so that should be biased as well..)
% % Hey Calvin, ACFs produced by your code are biased towards zero...
% % The reason for that is that the first k elements in cross_sum (variable of the sub-function) are always zero. Also, dimensions of cross_sum after the loop in lines 104-106 are always Nx1. In large sample the bias is small but in small samples it might be sensible.
% % Given that matlab is very bad at handling loops it is better to avoid them altogether if possible. I adjusted your code by removing the sub-function completely, "global" attributes for N and ybar (lines 46 and 48) and substituting loop in lines 52-54 by
% %
% % for i = 1:p
% % cross_sum=(y(i+1:N)-ybar)'*(y(1:N-i)-ybar);
% % yvar = (y-ybar)'*(y-ybar) ;
% % ta(i) = (cross_sum / yvar)*(N/(N-i)) ;
% % end
% %
% % Hope that helps everyone

% --------------------------
% USER INPUT CHECKS
% --------------------------

[n1, n2] = size(y);
if n2 ~=1
    error('Input series y must be an nx1 column vector')
end

[a1, a2] = size(nLags) ;
if ~((a1==1 && a2==1) && (nLags<n1))
    error('Input number of lags p must be a 1x1 scalar, and must be less than length of series y')
end



% -------------
% BEGIN CODE
% -------------

ta = zeros(nLags,1) ;
global N
N = max(size(y)) ;
global ybar
ybar = mean(y);

% Collect ACFs at each lag i
for i = 1:nLags
    ta(i) = acf_k(y,i) ;
end

lags    = 0:nLags;
ta      = [1; ta];

bounds(1,1) = (1.96)*(1/sqrt(N));
bounds(2,1) = (-1.96)*(1/sqrt(N));
numMA = 0;

% Plot ACF
if nargout == 0
    % Plot rejection region lines for test of individual autocorrelations
    % H_0: rho(tau) = 0 at alpha=.05
    
    lineHandles = stem(lags,ta,'filled','r-o');
    set(lineHandles(1),'MarkerSize',4)
    grid('on')
    xlabel('Lag')
    ylabel('Sample Autocorrelation')
    title('Sample Autocorrelation Function')
    hold('on')
    
    plot([numMA+0.5 numMA+0.5; nLags nLags],[bounds([1 1]) bounds([2 2])],'-b');
    
    plot([0 nLags],[0 0],'-k');
    hold('off')
    a = axis;
    axis([a(1:3) 1]);
else
    varargout = {ta, lags, bounds};
end

% line([0 p+.5], (1.96)*(1/sqrt(N))*ones(1,2))
% line([0 p+.5], (-1.96)*(1/sqrt(N))*ones(1,2))

% % Some figure properties
% line_hi = (1.96)*(1/sqrt(N))+.05;
% line_lo = -(1.96)*(1/sqrt(N))-.05;
% bar_hi = max(ta)+.05 ;
% bar_lo = -max(ta)-.05 ;
%
% if (abs(line_hi) > abs(bar_hi)) % if rejection lines might not appear on graph
%     axis([0 p+.60 line_lo line_hi])
% else
%     axis([0 p+.60 bar_lo bar_hi])
% end
% title({' ','Sample Autocorrelation Function',' '})
% xlabel('Lag')
% set(gca,'YTick',[-1:.20:1])
% % set number of lag labels shown
% if (p<28 && p>4)
%     set(gca,'XTick',floor(linspace(1,p,4)))
% elseif (p>=28)
%     set(gca,'XTick',floor(linspace(1,p,8)))
% end
% set(gca,'TickLength',[0 0])


% ---------------
% SUB FUNCTION
% ---------------
function ta2 = acf_k(y,k)
% ACF_K - Autocorrelation at Lag k
% acf(y,k)
%
% Inputs:
% y - series to compute acf for
% k - which lag to compute acf
%
global ybar
global N
cross_sum = zeros(N-k,1) ;

% Numerator, unscaled covariance
for i = (k+1):N
    cross_sum(i) = (y(i)-ybar)*(y(i-k)-ybar) ;
end

% Denominator, unscaled variance
yvar = (y-ybar)'*(y-ybar) ;

ta2 = sum(cross_sum) / yvar ;


