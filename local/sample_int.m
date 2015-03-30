function rand_int = sample_int(group0, group1)
% Sample from vectors containing indices
% (1) if the arguments are the same  -> assign a random element to each 
%     element(i) of _group0_ from _group0_ without considering the 
%     particular element(i) 
%
% (2) if the arguments are different -> assign a random element to each
%     element of _group0_ from _group1_
%
%USAGE:
% rand_int = SAMPLE_INT(group0, group1)
%
%INPUT:
% group0            sub-group of ensemble to make jump /vector of indices/
% group1            complementer sub-group of ensemble to get direction
%                   and step size for jump in case (1) /vector of indices/
%
%OUTPUT:
% rand_int         /vector/


if all(ismember(group1, group0)) %same arguments, sequential run
    Ngroup0         = numel(group0);
    rand_int        = randi(Ngroup0-1,1,Ngroup0);
    idx             = rand_int >= group0;
    rand_int(idx)   = rand_int(idx) + 1;   
else %different arguments, parallel run
    Ngroup0         = numel(group0);
    Ngroup1         = numel(group1);
    rand_int        = min(group1)-1 + randi(Ngroup1, 1, Ngroup0);
end

end