function [X0,mean0,C0] = initial_multiplicative(N,V,k,seed)
%INITIAL_MULTIPLICATIVE Initial data for the multiplicative-noise problem.
%   [X0,mean0,C0] = INITIAL_MULTIPLICATIVE(N,V,k,seed) returns the initial
%   ensemble, mean, and second moment.
%
%   Inputs
%   ------
%   N    : number of particles.
%   V    : orthonormal modal basis.
%   k    : number of active modes in the initial state.
%   seed : reserved input for deterministic initialization.
%
%   Outputs
%   -------
%   X0    : d-by-N initial ensemble.
%   mean0 : d-by-1 initial mean.
%   C0    : d-by-d non-centered second moment.

    s = linspace(0,1,size(V,1))';
    z = exp(-4*s(1:k));
    z = z/norm(z);

    mean0 = V(:,1:k)*z;
    X0 = mean0*ones(1,N);
    C0 = mean0*mean0';

end
