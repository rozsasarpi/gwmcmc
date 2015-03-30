function example_snowextremes
% Bayesian distribution fit to annual snow maxima
%
% Non-stationary generalized extreme value distribution (GEV)
%
% Xi ~ GEV(shape, scale, location_a*t + location_b)
%
% 't' is time, and the parameters has vague prior
%
% theta1 - shape
% theta2 - scale
% theta3 - location_a
% theta4 - location_b
%
% Datasource: CARPATCLIM Database © European Commission - JRC, 2013
%
% NOTE: This example depends on statistical toolbox!

close all
clc

% DATA
% snow water equivalent [mm] from 1962 to 2011; location ID: 300
x = [94, 93, 83, 126, 162, 106, 104, 79, 126, 125, 40, 64, 49, 21, 157,...
    19, 88, 125, 109, 61, 71, 17, 98, 95, 58, 95, 58, 68, 33, 39, 48, 41,...
    24, 22, 123, 28, 48, 116, 56, 39, 68, 62, 36, 102, 84, 52, 15, 18, 97, 36];

t       = 1:50;
n_obs   = numel(x);

% credible interval level
alpha_ci= 0.90;

%PRIOR
prior_shape         = @(xx) unicpdf(xx,-2,2);
prior_scale         = @(xx) unicpdf(xx,0,1000);
prior_location_a    = @(xx) normpdf(xx,0,1000);
prior_location_b    = @(xx) normpdf(xx,0,1000);

logprior = @(theta) log(prior_shape(theta(1)))+...
    log(prior_scale(theta(2)))+...
    log(prior_location_a(theta(3)))+...
    log(prior_location_b(theta(4)));

%LIKELIHOOD + PRIOR
logP_fun    = @(theta) logprior(theta) + loglike(theta);

nwalk = 1e2;    % for parallel run, increase the number of walks!
% set after an intial, exploring run (the parameters determine the upper bound
% of
link_init   = [
    (0.2*rand(1,nwalk)-0.1+1)*-0.1;
    (0.2*rand(1,nwalk)-0.1+1)*30;
    (0.2*rand(1,nwalk)-0.1+1)*-0.1;
    (0.2*rand(1,nwalk)-0.1+1)*80];

Options.parallel = 0;   %compare the parallel and sequential analyses!
%MCMC hammer
tic
%--------------------------------------------------------------------------
ensembles   = gwmcmc_par(link_init, logP_fun, 1e5, Options);
% ensembles   = gwmcmc0(link_init, {logP_fun}, 1e5);
%--------------------------------------------------------------------------
toc

% burn-in
ensembles(:,:,1:20) = [];

Plot.autocorr = 1;
% Plot.triangle = 1;
num_table = gwmcmc_diag(ensembles, Plot);
num_table.Properties.RowNames = {'shape', 'scale', 'location_a', 'location_b'};
disp(num_table)

% %PLOT
RP      = [50, logspace(log10(1.01),log10(1900),10)];

RP      = unique(RP);
P       = 1 - 1./RP;
empiF   = ((1:(length(x)))-2/5)/(length(x)+1/5);
% no general agreement on how to construct empirical distribution, affects only the plot!

squeezed_ens = ensembles(:,:)';
% year for plot (1962+tt)
tt = 25;
shape_par   = squeezed_ens(:,1);
scale_par   = squeezed_ens(:,2);
loc_par     = squeezed_ens(:,3)*tt + squeezed_ens(:,4);
q_ci        = nan(3, numel(RP));
for jj = 1:numel(RP)
    q = gevinv(repmat(P(jj),size(shape_par)),shape_par, scale_par,loc_par);
    q_ci(:,jj) = quantile(q, [(1-alpha_ci)/2, 0.5, (1+alpha_ci)/2]);    
end

tRP = -log(-log(1-1./RP));
figure('Position',[200, 200, 600, 400]);

plot(-log(-log(empiF)), sort(x),'o',...
    'Color','black',...
    'MarkerSize',3);
hold on
plot(tRP, q_ci(2,:), 'Color', [1 0 0], 'LineWidth', 1.25);
plot(tRP, q_ci(1,:), 'Color', [1 0 0], 'LineStyle', '--', 'LineWidth', 0.75);
plot(tRP, q_ci(3,:), 'Color', [1 0 0], 'LineStyle', '--', 'LineWidth', 0.75);

xlabel('Return period [year]')
ylabel('SWE [mm]')

xmin = min(-log(-log(empiF)))- ( max(-log(-log(empiF))) - min(-log(-log(empiF))) )*0.05;%, -2.3212*1.1);
xmax = max(max(-log(-log(empiF)))*1.7, 6.9073*1.1);
ymin = min(x) - 0.05*(max(x)-min(x));
ymax = max(x) + 0.6*(max(x)-min(x));
axis([xmin xmax ymin ymax])

xtick = [0.1 1 10 100 1000]; % WARNING all should be power of 10!
xtick_tmp = -log(-log(1-1./xtick(end-1:end)));
xtick_diff = diff(xtick_tmp);
xtick_transf = xtick_tmp(2) - xtick_diff*(0:length(xtick)-1);
xtick_transf = fliplr(xtick_transf);
set(gca,'XTick',xtick_transf);
set(gca,'XTickLabel',xtick);

title(['\bf Non-stationary GEV fitted to annual snow maxima, CI:',num2str(alpha_ci),', tt = ', num2str(tt)])
grid on
%==========================================================================
% NESTED FUNCTIONS
%==========================================================================
    function LL = loglike(theta)
        shape       = theta(1);
        scale       = theta(2);
        location_a  = theta(3);
        location_b  = theta(4);
        
        shape_      = repmat(shape,1,n_obs);
        scale_      = repmat(scale,1,n_obs);
        location_   = location_a*t + location_b;
        
        LL = sum(log(gevpdf(x, shape_, scale_, location_)));        
    end

    function f = unicpdf(xx, a, b)
        f          = nan(size(xx));
        idx        = xx<a | xx>b;
        f(idx)     = 0;
        f(~idx)    = 1/(b-a);
    end

end