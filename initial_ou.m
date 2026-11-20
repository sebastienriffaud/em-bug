function [X0,mean0,C0] = initial_ou(d,N,seed)
%INITIAL_OU Initial data for the OU problem.
%   [X0,mean0,C0] = INITIAL_OU(d,N,seed) returns the initial
%   ensemble, mean, and second moment.
%
%   Inputs
%   ------
%   d    : state dimension.
%   N    : number of particles.
%   seed : reserved input for deterministic initialization.
%
%   Outputs
%   -------
%   X0    : d-by-N initial ensemble.
%   mean0 : d-by-1 initial mean.
%   C0    : d-by-d initial non-centered second moment.

    mean0 = zeros(d,1);
    mean0(1) = 1;
    X0 = mean0*ones(1,N);
    C0 = mean0*mean0';

end
