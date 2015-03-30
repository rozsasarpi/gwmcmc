function logP = get_logP(link, Options)
% Evaluates the logP_fun function
%
%USAGE:
% logP = GET_LOGP(link, Global)
%
%INPUT:
% link          _link_ positions where the _logP_ is needed /MxN matrix/
%               M (dim) is the number of inferred parameters
% Global        /structure/
%  .logP_fun    /function handler/
%
%OUTPUT:
% logP      /1xN matrix/

parallel    = Options.parallel;
logP_fun    = Options.logP_fun;
Nlink       = size(link,2);

logP        = nan(1,size(link,2));

%loop over the links in the walkers
switch parallel
    case 1
        parfor ii = 1:Nlink
            theta       = link(:,ii);
            logP(ii)    = logP_fun(theta);
        end
        
    case 0
        % this for loop is used only during the initialization
        % during the MCMC simulation the loop is outside of this function
        % >> to reflect the moves of the previous walkers in the current
        for ii = 1:Nlink
            theta       = link(:,ii);
            logP(ii)    = logP_fun(theta);
        end
end

end