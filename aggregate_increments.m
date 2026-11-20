function dWcoarse = aggregate_increments(dWfine,factor)
%AGGREGATE_INCREMENTS Sum fine-grid Brownian increments over coarse steps.
%   dWcoarse = AGGREGATE_INCREMENTS(dWfine,factor) groups consecutive
%   increments in the third dimension and sums each group.
%
%   Inputs
%   ------
%   dWfine : q-by-N-by-Ntfine array of Brownian increments.
%   factor : positive integer, number of fine steps per coarse step.
%
%   Output
%   ------
%   dWcoarse : q-by-N-by-(Ntfine/factor) array.
%
%   Ntfine must be divisible by factor.

    [q,N,Ntfine] = size(dWfine);
    % Brownian increments add across consecutive time intervals.
    Ntcoarse = Ntfine/factor;
    dWcoarse = zeros(q,N,Ntcoarse);

    for k = 1:Ntcoarse
        indices = (k-1)*factor+(1:factor);
        dWcoarse(:,:,k) = sum(dWfine(:,:,indices),3);
    end

end
