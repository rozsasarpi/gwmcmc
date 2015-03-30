function [ensembles, logP, Options] = gwmcmc_ini(link_init, logP_fun, mc_count, Options)
% Initialization of gwmcmc sampler.
%
% GWMCMC_INI

% set default values if option fields are not defined
Option_fields   = {'skip', 'stepsize', 'silent',    'parallel', 'paropti'};
Option_defaults = [10,      2.5,        0,          0,          0];         %WARNING! order!

if nargin < 3
    error('GWMCMC requires at least 3 inputs.')
end

dim         = size(link_init,1);

if size(link_init,2) == 1
    link_init=bsxfun(@plus,link_init,randn(dim,dim*5));
end

nwalkers = size(link_init,2);

for ii = 1:numel(Option_fields)
    if (nargin < 4) || ~isfield(Options, Option_fields{ii})
        Options.(Option_fields{ii}) = Option_defaults(ii);
    end
end

skip        = Options.skip;
stepsize    = Options.stepsize;
silent      = Options.silent;
parallel    = Options.parallel;
paropti     = Options.paropti;

% _Options_ structure variable, carries the fix variables (for more concise argument passing)
Options.logP_fun = logP_fun;
Options.dim      = dim;
Options.nwalkers = nwalkers;

if paropti == 1 && parallel == 0
    warning('MATLAB:input',...
        'Paropti is true but parallel is false (or not specified)!\n >> Sequential run!')
end

if paropti == 1 && parallel == 1
    warning('MATLAB:Info',...
        'The provided ''link_init'' is overwritten to optimize the parallel run!\n See the help for details!')
    % here comes the function which overwrites
end

if size(link_init,1)*2 > size(link_init,2)
    warning('MATLAB:input',...
        'Check ''link_init'' dimensions!\n It is recommended that there be at least twice as many walkers in the ensemble as there are model dimensions. ')
end

Nkeep       = ceil(mc_count/skip/nwalkers);  %number of samples drawn from each walker
mc_count    = max((Nkeep-1)*skip+1,2);

ensembles   = nan(dim,nwalkers,Nkeep);        %pre-allocate output matrix

ensembles(:,:,1) = link_init;

if iscell(logP_fun)
    error('The cell type logP_fun defintion is no longer supported, see gwmcmc help!')
end

%calculate logP state initial pos of walkers
logP        = nan(1,nwalkers,Nkeep);
logP(1,:,1) = get_logP(link_init, Options);

% WARNING!
if ~all(isfinite(logP(:,:,1)))
    error('Starting points for all walkers must have finite logP!')
end

Options.mc_count = mc_count;

end