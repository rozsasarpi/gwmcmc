% Dugongs: a nonlinear growth curve
%
% This example illustrates a simple nonlinear, Bayesian regression.
% Example adopted from: http://www.openbugs.net/Examples/Dugongs.html
% original: Carlin & Gelfand.1991. An iterative Monte Carlo method for nonconjugate Bayesian analysis
%
% Y_i ~ Normal(mu_i, tau)
%
% mu_i = alpha - beta*gamma^x_i;    alpha, beta > 0; 0 < gamma < 1
%
% the parameters has vague prior
% theta1 - alpha
% theta2 - beta
% theta3 - gamma
% theta4 - sigma    %instead of precision standard deviation is used directly!
%
% theta5 - U3       %not required to define the model!, will be derived from _gamma_
%
% NOTE: This example is still too simple to really take advantage of
%       parallel computing, although the runtimes are comparable (tested on 4-core processor)

function example_dugongs

% clear all
close all
clc
rng(125)

%DATA
x = [1.0,  1.5,  1.5,  1.5, 2.5,   4.0,  5.0,  5.0,  7.0, 8.0,  8.5,  9.0,  9.5, 9.5,...
    10.0, 12.0, 12.0, 13.0, 13.0, 14.5, 15.5, 15.5, 16.5, 17.0, 22.5, 29.0, 31.5];
y = [1.80, 1.85, 1.87, 1.77, 2.02, 2.27, 2.15, 2.26, 2.47, 2.19, 2.26, 2.40, 2.39, 2.41,...
    2.50, 2.32, 2.32, 2.43, 2.47, 2.56, 2.65, 2.47, 2.64, 2.56, 2.70, 2.72, 2.57];

n_obs = numel(x);

%PRIOR
prior_alpha = @(xx) trunc_normpdf(xx,0,1000);
prior_beta  = @(xx) trunc_normpdf(xx,0,1000);
prior_gamma = @(xx) unicpdf(xx,0.5,1.0);
prior_sigma = @(xx) unicpdf(xx,0,1000);
% prior_tau   = @(xx) gampdf(xx,1e-3,1e3);

logprior = @(theta) log(prior_alpha(theta(1)))+...
    log(prior_beta(theta(2)))+...
    log(prior_gamma(theta(3)))+...
    log(prior_sigma(theta(4)));

%LIKELIHOOD + PRIOR
logP_fun    = @(theta) logprior(theta) + loglike(theta);

nwalk = 1e2;    % for parallel run, increase the number of walks!
% set after an intial, exploring run
link_init   = [(0.1*rand(1,nwalk)+0.95)*2.65;
    (0.1*rand(1,nwalk)+0.95)*0.97;
    (0.1*rand(1,nwalk)+0.95)*0.86;
    (1*rand(1,nwalk)+0.50)*0.1];

% for exploring
% link_init   = [100*rand(1,nwalk);
%                100*rand(1,nwalk);
%                0.5*rand(1,nwalk)+0.5;
%                100*rand(1,nwalk)];

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

% add U3 variable to _ensembles_ to get diagnostic of it
% U3 = log(gamma/(1+gamma))
npar = size(ensembles,1);
ensembles(npar+1,:,:) = log(ensembles(3,:,:)./(1-ensembles(3,:,:)) );

gwmcmc_diag(ensembles)

% PLOT the fitted model
%--------------------------------------------------------------------------
% squeeze
squeezed_ens = ensembles(:,:)';

dage = linspace(min(x),max(x),20);
ndage = numel(dage);
dlength = nan(3,ndage);

for jj = 1:ndage
    dlength_sample = squeezed_ens(:,1) - squeezed_ens(:,2).*squeezed_ens(:,3).^dage(jj);
    dlength(:,jj)  = quan_tile(dlength_sample, [0.025, 0.50, 0.975]);
end

figure
% !the credible interval represents only the uncertainty in the mean(mu) estimate!!
plot(dage, dlength([1,3],:),'b--')
hold on
plot(dage, dlength(2,:), 'Color', 'b', 'LineWidth', 1.5)
plot(x,y, 'or')

img = imread('dugong.png');
img = imrotate(img,180);
image([20,33],[1.7,2.2], img) 
xlabel('Age [year]')
ylabel('Length [m]')
title('\bf Dugong growth curve')
grid on

%==========================================================================
% NESTED FUNCTIONS
%==========================================================================
    function LL = loglike(theta)
        alpha   = theta(1);
        beta    = theta(2);
        gamma   = theta(3);
        sigma   = theta(4);
        
        mu      = alpha - beta*gamma.^x;
        
        LL = 0;
        for ii = 1:n_obs
            LL = log(normpdf(y(ii), mu(ii), sigma)) + LL;
%             if ~isfinite(LL)  %WARNING!
%                 LL = -1e100;
%                 break
%             end
        end
        
    end

    function f = trunc_normpdf(x,mu,sigma) %[0,Inf]
        f          = nan(size(x));
        idx        = x<0;
        f(idx)     = 0;
%         f(~idx)    = 2*normpdf(x(~idx),mu,sigma);
        f(~idx)    = 1/( sqrt(2*pi) * sigma ) * exp( -1/2 * ((x(~idx)-mu)/sigma).^2 );
    end

    function f = unicpdf(xx, a, b)
        f          = nan(size(xx));
        idx        = xx<a | xx>b;
        f(idx)     = 0;
        f(~idx)    = 1/(b-a);
    end

end