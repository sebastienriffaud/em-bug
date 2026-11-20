function [drift,noise,G,H] = ou(~,X,C,A,B,dW)
%OU Ornstein--Uhlenbeck problem operators.
%   [drift,noise,G,H] = OU(t,X,C,A,B,dW) returns the drift, noise, and
%   moment operators.
%
%   Inputs
%   ------
%   t  : time (unused; this problem is autonomous).
%   X  : d-by-N ensemble.
%   C  : d-by-d non-centered second moment.
%   A  : d-by-d drift matrix.
%   B  : d-by-m diffusion matrix.
%   dW : m-by-N Brownian increments.
%
%   Outputs
%   -------
%   drift, noise : d-by-N Euler--Maruyama updates.
%   G, H         : d-by-d moment-update matrices.

    drift = A*X;
    noise = B*dW;

    if nargout > 2
        Q = B*B';
        G = A*C+C*A'+Q;
        H = A*C*A';
    end

end
