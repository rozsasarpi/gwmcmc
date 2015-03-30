function [ensembles, logP] = gwmcmc_par(link_init, logP_fun, mc_count, Options)
%% Affine invariant ensemble MCMC sampler. "The MCMC hammer"
%
% Ensemble Markov Chain Monte Carlo sampling of posterior distribution
% Only for continuous functions!
% (A variant of the Goodman and Weare affine invariant ensemble sampler).
%
%USAGE:
%  [ensembles, logP] = GWMCMC(link_init, logPf_uns, mc_count, options);
%
%INPUTS:
%  link_init:   an MxW matrix of initial values for each of the walkers in the
%               ensemble. (M:number of model params. W: number of walkers). W
%               should be atleast 2xM. (see e.g. mvnrnd) /MxW matrix/.
%
%  logP_fun:    loglikelihood function with input as the vector of parameters,
%               to pass prior: logP_fun = @(theta) loglike(theta) + logprior(theta)
%               /function handle/.
%
%  mc_count:    What is the desired total number of monte carlo proposals.
%               This is the total number, -NOT the number per chain /integer/.
%
%  Options:     optional! control /structure/
%   .stepsize:  unit-less stepsize /integer/
%               [default=2.5].
%   .skip:      Thin all the chains by only storing every N'th step /integer/
%               [default=10].
%   .silent:    Run in silent mode, without progress report /logical, integer/
%               [default=0].
%   .parallel:  Run the walker updating in parallel for each simulation step (chain link)!
%               (the ensemble is split into two subgroups to maintain the detailed balance)
%               This option is only recommended if _logP_fun_ is costly because
%               the overhead in distributing and assembling jobs between the
%               workers (see tips for more detailes) /logical, integer/
%               [default=0].
%   .paropti:   Optimize the number of walkers (partitioning the jobs) if _parallel_
%               option is selected. If this set to true then the provided _link_init_
%               is overwritten (it is used to construct the 
%               optimized _link_init_); /logical, integer/
%               [default=0]. !!not yet implemented!!
%
%OUTPUTS:
%   ensembles:  A MxWxT matrix with the thinned markov chains (with T samples
%               per walker). T=~mccount/skip/W.
%   logP:       A 1xWxT matrix of log probabilities for each model in the
%               ensembles.
%
%CUSTOM FUNCTIONS: (only the directly called)
%  propose_stretch.m     Propose new position for one subgroup of ensemble
%  get_logP.m            Evaluates the _logP_fun_ function (parallel)

% Note the cascaded evaluation is no longer supported. TODO: add again
%
% TIPS:(1)  If you aim to analyze the entire set of ensemble members as a single
%           sample from the distribution then you may collapse output models-matrix
%           thus: models=models(:,:); This will reshape the MxWxT matrix into a
%           Mx(W*T)-matrix while preserving the order.
%      (2)  Parallel run: the chains must evolve sequentially but the links of ensembles
%           can be evaluated parallel. To maintain the detailed balance,
%           the links are splitted into two subgroups. Each group is evaluated
%           in parallel but the groups are sequentially following each other.
%           The second group of links are moved based on the move of the first
%           group. If the number of chain links is large && few walkers
%           are applied && _logP_fun_ is cheap to evaluate, the parallel calculation
%           will likely take longer then the sequential run!
%           Thus the parallel run is only beneficial if the overhead is counterbalanced
%           by the gain during the parallel calculation (reduce chain links &&
%           increase ensemble number=walkers).
%           When more processors(threads) are available the parallel option can be
%           faster for smaller problems as well.
%           Furthermore, it is recommended to use _.paropti_ to roughly optimze
%           the number of chain links and ensembles.
%
% EXAMPLE:
%   %Here we sample a multivariate normal distribution.
%   %define problem:
%   mu = [5;-3;6];
%   C  = [.5 -.4 0;-.4 .5 0; 0 0 1];
%   iC = pinv(C);
%   logP_fun = @(m)-0.5*sum((m-mu)'*iC*(m-mu))
%
%   %make a set of starting points for all the ensemble of walkers
%   link_init = randn(length(mu),length(mu)*2);
%
%   %Apply the MCMC hammer
%   [ensembles,logP] = GWMCMC(link_init,logP_fun,1e5);
%   burnin = floor(size(ensembles,3)/5);
%   ensembles(:,:,1:burnin) = []; %remove 20% as burn-in
%   ensembles = ensembles(:,:)';           %reshape matrix to collapse the ensemble member dimension
%   scatter(ensembles(:,1),ensembles(:,2))
%   prctile(ensembles,[5 50 95])
%
%
% References:
% Goodman & Weare (2010), Ensemble Samplers With Affine Invariance,
%    Comm. App. Math. Comp. Sci., Vol. 5, No. 1, 65–80
% Foreman-Mackey, Hogg, Lang, Goodman (2013), emcee: The MCMC Hammer,
%    arXiv:1202.3665
%
% See also
% gwmcmc_demo, gwmcmc_diag
%
% -Aslak Grinsted 2015

%--------------------------------------------------------------------------
%% INPUT CHECK AND INITIALIZATION
%--------------------------------------------------------------------------
if exist('Options','var')
    [ensembles, logP, Options] = gwmcmc_ini(link_init, logP_fun, mc_count, Options);
else
    [ensembles, logP, Options] = gwmcmc_ini(link_init, logP_fun, mc_count);
end

% basic control parameters
mc_count        = Options.mc_count;
silent          = Options.silent;
skip            = Options.skip;
nwalkers        = Options.nwalkers;
parallel        = Options.parallel;

% set the first links
current_link    = ensembles(:,:,1);
current_logP    = logP(:,:,1);

% progress bar
if silent == 0
    ctime       = cputime;
    starttime   = cputime;
    progress(0,current_link,0)
else
    %do nothing
end

%split the ensemble into two subgroups (for parallel computing)
split_pos   = floor(nwalkers/2);
group0      = 1:split_pos;
group1      = (split_pos + 1):nwalkers;
group_tot   = 1:nwalkers;

%--------------------------------------------------------------------------
%% MCMC SIMULATION
%--------------------------------------------------------------------------
% loop over the links
reject = 0;
for ii = 2:mc_count
    
    switch parallel
        case 1
            current_logP0 = current_logP(group0);
            current_logP1 = current_logP(group1);
            
            % move group0 based on group1
            [new_link0, new_logP0, accept0] = propose_stretch(group0, group1, current_link, current_logP0, Options);
            
            %update the link pos and logP of group0
            current_link(:,group0) = new_link0;
            current_logP(:,group0) = new_logP0;
            
            % move group1 based on the updated positions of group0
            [new_link1, new_logP1, accept1] = propose_stretch(group1, group0, current_link, current_logP1, Options);
            
            %update the link pos and logP of group1
            current_link(:,group1) = new_link1;
            current_logP(:,group1) = new_logP1;
            
            accept = [accept0, accept1];
        case 0
            [new_link, new_logP, accept] = propose_stretch(group_tot, group_tot, current_link, current_logP, Options);
       
            %update the link pos and logP
            current_link = new_link;
            current_logP = new_logP;
    end
    
    %count rejected proposals
    reject = reject + sum(not(accept));
    
    % thinning
    if mod(ii-1,skip) == 0
        row                 = ceil(ii/skip);
        ensembles(:,:,row)  = current_link;
        logP(1,:,row)       = current_logP;
    end
    
    %progress bar
    if silent == 0
        progress(ii/mc_count, mean(current_link,2), reject/(ii*nwalkers))
    else
        %do nothing
    end
    
end

%progress bar
if silent == 0
    progress(1, mean(current_link,2), reject/(ii*nwalkers));
else
    %do nothing
end


% TODO: make standard diagnostics to give warnings...
% TODO: cut away initial drift.(?)
% TODO: make some diagnostic plots if nargout==0;

% TODO: initialization and preprocess in separate m-file
%--------------------------------------------------------------------------
%% NESTED FUNCTIONS
%--------------------------------------------------------------------------
function progress(pct,curm,rejectpct)
% Progress report
persistent lastNchar lasttime starttime
if isempty(lastNchar) || pct == 0
    lastNchar = 0;lasttime = cputime-10;starttime = cputime;fprintf('\n')
    pct = 1e-16;
end
if pct == 1
    fprintf('%s',repmat(char(8),1,lastNchar));lastNchar = 0;
    return
end
if (cputime-lasttime>0.1)
    
    ETA = datestr((cputime-starttime)*(1-pct)/(pct*60*60*24),13);
    progressmsg = [uint8((1:40) <= (pct*40)).*'#' ''];
    curmtxt = sprintf('% 9.3g\n',curm(1:min(end,20),1));
    %curmtxt = mat2str(curm);
    progressmsg = sprintf('GWMCMC %5.1f%% [%s] %s\n%3.0f%% rejected\n%s\n',pct*100,progressmsg,ETA,rejectpct*100,curmtxt);
    
    fprintf('%s%s',repmat(char(8),1,lastNchar),progressmsg);
    drawnow;lasttime = cputime;
    lastNchar = length(progressmsg);
end

% Acknowledgements: I became aware of the algorithm via a student report
% which was using emcee for python. I read the paper and judged that this
% must be excellent, and made my own implementation for matlab. It is.

% Inspired by the emcee python implementation
% http://dan.iel.fm/emcee/current/