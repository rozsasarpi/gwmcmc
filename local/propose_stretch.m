function [new_link0, new_logP0, accept0] = propose_stretch(group0, group1, current_link, current_logP0, Options)
% Propose new position for one subgroup of ensemble
%
%USAGE:
%[new_link0, new_logP0, accept0] = PROPOSE_STRETCH(group0, group1, current_link, current_logP0, Global)
%
%INPUT:
% group0            sub-group of ensemble to make jump /vector of indices/
% group1            complementer sub-group of ensemble to get direction
%                   and step size for jump /vector of indices/
%
% current_link
% current_logP0
% Options
%
%OUTPUT:
% new_link0         /vector/
% new_logP0         /vector/
% accept0           /boolen/
%
%CUSTOM FUNCTIONS: (only the directly called)
% get_logP.m        Evaluates the _logP_fun_ function
% sample_int.m      Sample from vectors containing indices
%
%NOTE:  Works with sequential running as well, in that case group0 has one
%       element only.


stepsize            = Options.stepsize;
dim                 = Options.dim;

Ngroup0             = numel(group0);
rand_int            = sample_int(group0, group1);
current_link0       = current_link(:,group0);
rand_link1          = current_link(:,rand_int);

zz                  = ((stepsize - 1)*rand(1,Ngroup0) + 1).^2/stepsize;
proposed_link0      = rand_link1 + bsxfun(@times, (current_link0 - rand_link1), zz);
proposed_logP0      = get_logP(proposed_link0, Options);

q                   = (dim-1)*log(zz) + proposed_logP0 - current_logP0;

accept0             = log(rand(1,Ngroup0)) < q;

new_link0           = current_link0;
new_link0(:,accept0)  = proposed_link0(:,accept0);

new_logP0           = current_logP0;
new_logP0(accept0)  = proposed_logP0(accept0);

end

%TODO: input and output could be more concise, but parfor doesnt like structure variables..
%
% Snippet to break the comparison into two phases, to ulitilize if the cheap
% logPfuns is sufficient to reject, rarely happens, is it worth to keep this system?
%         % supremum on _q_ (try if the step can be rejected only by using
%         %the cheaper _logPfuns_)
%         sup_q = (M-1)*log(zz) + proposed_logP1 - current_logP(1,w_idx);
%         if logrand < sup_q
%             % calculate the full _q_, add the more expensive _logPfuns_ as well
%             proposedlogP2 = logP_funs{2}(proposed_link);
%             if logrand < sup_q + proposedlogP2
%                 proposed_logP    = proposed_logP1 + proposedlogP2;
%             else
%                 reject          = reject + 1;
%                 acceptstep      = false;
%             end
%         else
%             reject          = reject + 1;
%             acceptstep      = false;
%         end