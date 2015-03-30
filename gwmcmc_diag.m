function num_table = gwmcmc_diag(ensembles, Plot)
%% Diagnostic plot for gwmcmc results
%
%USAGE:
% GWMCMC_DIAG(ensembles, Plot)
%
%INPUT:
% ensembles     3 dimensional matrix, output of gwmcmc, Npar x Nwalk x Nlink 
%OPTIONAL:
% Plot          specifying the desired figures /structure/
%  .autocorr    autocorrelation plots per walker /boolen, integer, 0 or 1/
%  .chain       trace plot of chains per parameter /boolen, integer, 0 or 1/
%  .triangle    triangle plot with 1D and 2D marginals of the parameters,
%               correlations and credible intervals /boolen, integer, 0 or 1/
%  .hist        1D marginal of the parameters,
%               median, mean and credible intervals /boolen, integer, 0 or 1/
%
%OUTPUT:
% num_table     table with numerical summary
%
%CUSTOM FUNCTIONS:
% diag_*.m      (separate m-file for every diagnostic plots, with name 
%               matching that of used in _Plot_
% quan_tile.m
%
% See also
% gwmcmc

%TODO: kernel density plot (kde by Botev)
%TODO: check if _Plot_ fields are mistyped

%--------------------------------------------------------------------------
% Input check and initialization
%--------------------------------------------------------------------------
if nargin < 1
    error('''gwmcmc_diag'' requires at least 1 input.')
end

% set default values if option fields are not defined
Plot_fields   = {'autocorr', 'chain', 'triangle', 'hist'};
% Plot_defaults = [1,           1,       1,             1];         %WARNING! order!
Plot_defaults = [0,           0,       1,             0];

% set defaults if not specified
for ii = 1:numel(Plot_fields)
    if (nargin < 2) || ~isfield(Plot, Plot_fields{ii})
        Plot.(Plot_fields{ii}) = Plot_defaults(ii);
    end
end

if numel(size(ensembles)) < 3
    error('Ensembles should be a 3-dimensional matrix, output of gwmcmc!')
end

npar = size(ensembles,1);

%--------------------------------------------------------------------------
% Plot diagnostic figures
%--------------------------------------------------------------------------
for ii = 1:numel(Plot_fields)
    if Plot.(Plot_fields{ii}) == 1
        feval(['diag_',Plot_fields{ii}], ensembles)
    end
end

%--------------------------------------------------------------------------
% Numeric summary
%--------------------------------------------------------------------------
% squeeze
squeezed_ens = ensembles(:,:)';

mean_theta      = mean(squeezed_ens)';
median_theta    = median(squeezed_ens)';
sd_theta        = std(squeezed_ens)';
q_theta         = quan_tile(squeezed_ens,[0.025, 0.975])';
sample_theta    = repmat(size(squeezed_ens,1), npar, 1);
MC_sd_error     = sd_theta./sqrt(sample_theta); %?! not sure (effect of thinning etc.)

num_table       = table(mean_theta, sd_theta, MC_sd_error, q_theta(:,1), median_theta, q_theta(:,2), sample_theta);
% keyboard
num_table.Properties.VariableNames = {'mean', 'sd', 'MC_sd_error', 'q_025', 'median_q_500', 'q_975', 'sample'};

row_names = cell(npar,1);
for ii = 1:npar
    row_names{ii} = ['theta_',num2str(ii)];
end
num_table.Properties.RowNames = row_names;

end