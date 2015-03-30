% Rats: a normal hierarchical model
%
% Example adopted from: http://www.openbugs.net/Examples/Rats.html
%
% Note: The notation follows that of the OpenBUGS' example.
%       Instead of the precision parameters, standard deviation is used.
%       Vague, uniform prior is adopted for the standard deviations.
%
% theta1 - alpha_c
% theta2 - alpha_sigma
% theta3 - beta_c
% theta4 - beta_sigma
% theta5 - sigma_c
%
% theta6:35  - alpha..
% theta36:65 - beta..
%
% Rat weight - day relationship
% Each rat has its own linear relationship, however these are not independent
% but connected (constrained) by common hyperprior distributions.
%
% NOT WORKING PROPERLY YET!


function example_rats

close all
% clear all
clc

%DATA
Y = [151   145   147   155   135   159   141   159   177   134   160   143   154   171   163   160   142   156   157   152   154   139   146   157   132   160   169   157   137   153
    199   199   214   200   188   210   189   201   236   182   208   188   200   221   216   207   187   203   212   203   205   190   191   211   185   207   216   205   180   200
    246   249   263   237   230   252   231   248   285   220   261   220   244   270   242   248   234   243   259   246   253   225   229   250   237   257   261   248   219   244
    283   293   312   272   280   298   275   297   350   260   313   273   289   326   281   288   280   283   307   286   298   267   272   285   286   303   295   289   258   286
    320   354   328   297   323   331   305   338   376   296   352   314   325   358   312   324   316   317   336   321   334   302   302   323   331   345   333   316   291   324];

Y       = Y';
% Y       = Y(1:10,:);
x       = [8, 15, 22, 29, 36];
x_bar   = 22;

[nrat, nday] = size(Y);
npar    = 5+2*nrat;

%PRIOR (~hyperprior)
prior_alpha_c       = @(xx) normpdf(xx,0,1000);
prior_alpha_sigma   = @(xx) unicpdf(xx,0,10000);
prior_beta_c        = @(xx) normpdf(xx,0,1000);
prior_beta_sigma    = @(xx) unicpdf(xx,0,10000);
prior_sigma_c       = @(xx) unicpdf(xx,0,1000);

logprior = @(theta) log(prior_alpha_c(theta(1)))+...
    log(prior_alpha_sigma(theta(2)))+...
    log(prior_beta_c(theta(3)))+...
    log(prior_beta_sigma(theta(4)))+...
    log(prior_sigma_c(theta(5)));

%LIKELIHOOD + PRIOR
logP_fun    = @(theta) logprior(theta) + loglike(theta);

% MLE estimate
% options = optimset('MaxFunEvals', 1e4, ...
%             'MaxIter', 1e4, ...
%             'TolFun', 1e-6, ...
%             'TolX', 1e-6);
% theta_ini = [215.0951    8.8094    7.3908   12.0670    3.959];
% theta_ini = [theta_ini, ones(1,nrat)*239.8, ones(1,nrat)*6.01];
% [theta_mle, maxLL] = fminsearch(@(theta) -loglike(theta), theta_ini, options);
% % [theta_mle, maxLL] = ga(@(theta) -loglike(theta), numel(theta_ini));
% 
% figure
% plot(x,Y,'--b')
% hold on
% a = theta_mle(1);
% b = theta_mle(3);
% a1 = theta_mle(6:(6+nrat-1));
% b1 = theta_mle((6+nrat):(6+2*nrat-1));
% for u = 1:numel(a1)
%     plot(x,a1(u)+b1(u)*(x-x_bar),'g')
% end
% plot(x,a+b*(x-x_bar),'Color','r','LineWidth', 2)

nwalk = 1e2;
% set after an intial run
link_init   = nan(npar,nwalk);
link_init(1:5,:)   = [
    (0.2*rand(1,nwalk)-0.1+1)*242;       %alpha_c
    (0.2*rand(1,nwalk)-0.1+1)*21;            %alpha_sigma
    (0.2*rand(1,nwalk)-0.1+1)*6;               %beta_c
    (0.2*rand(1,nwalk)-0.1+1)*2;            %beta_sigma
    (0.2*rand(1,nwalk)-0.1+1)*19];           %sigma_c

link_init(6:(6+nrat-1),:)           = (0.2*rand(nrat,nwalk)-0.1+1)*210;     %alpha
link_init((6+nrat):(6+2*nrat-1),:)  = (0.2*rand(nrat,nwalk)-0.1+1)*6;     %beta

% Options.stepsize = 2;
Options.parallel = 0;   %compare the parallel and sequential analysis
%MCMC hammer
tic
%--------------------------------------------------------------------------
[ensembles, logP]   = gwmcmc_par(link_init, logP_fun, 1e5, Options);
% ensembles   = gwmcmc0(link_init, {logP_fun}, 1e5);
%--------------------------------------------------------------------------
toc

% burn-in
ensembles(:,:,1:20) = [];

% gwmcmc_diag(ensembles)

% add alpha_0, intercept
% ensembles(npar+1,:,:) = ensembles(1,:,:) - ensembles(3,:,:)*x_bar;
% keep only the intercept, slope and standard deviation for diagnostics
% ensembles = ensembles([npar+1,3,7],:,:);
ensembles = ensembles(1:7,:,:);

num_table = gwmcmc_diag(ensembles);
% num_table.Properties.RowNames = {'alpha_0', 'beta_c', 'sigma_c'};
disp(num_table)

%PLOT
squeezed_ens = ensembles(:,:)';
figure
plot(x,Y,'--b')
hold on
a = mean(squeezed_ens(:,1));
b = mean(squeezed_ens(:,3));
plot(x,a+b*(x-x_bar),'Color','r','LineWidth',2)
grid on

disp('completed')
%==========================================================================
% NESTED FUNCTIONS
%==========================================================================
    function LL = loglike(theta)
        alpha_c         = theta(1);
        alpha_sigma     = theta(2);
        beta_c          = theta(3);
        beta_sigma      = theta(4);
        sigma_c         = theta(5);
        
        % for each rat
        alpha           = theta(6:(6+nrat-1));
        beta            = theta((6+nrat):(6+2*nrat-1));
                
        sigma   = repmat(sigma_c,1,nday);
        
        % this is like the first level prior
        LL = sum(log(normpdf(alpha, alpha_c, alpha_sigma)))+...
            sum(log(normpdf(beta, beta_c, beta_sigma)));

        %loop over the rats
        for ii = 1:nrat
            mu = alpha(ii) + beta(ii)*(x - x_bar);
            LL = sum(log(normpdf(Y(ii,:),mu,sigma))) + LL;
            %             if ~isfinite(LL)  %WARNING!
            %                 LL = -1e100;
            %                 break
            %             end
        end
    end

    function f = unicpdf(xx, a, b)
        f          = nan(size(xx));
        idx        = xx<a | xx>b;
        f(idx)     = 0;
        f(~idx)    = 1/(b-a);
    end

end