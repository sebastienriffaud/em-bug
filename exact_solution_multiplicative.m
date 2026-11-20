function X = exact_solution_multiplicative(X0,V,lambda,mu,dt,dW)
%EXACT_SOLUTION_MULTIPLICATIVE Exact ensemble for the multiplicative problem.
%   X = EXACT_SOLUTION_MULTIPLICATIVE(X0,V,lambda,mu,dt,dW) returns the
%   exact ensemble.
%
%   Inputs
%   ------
%   X0        : d-by-N initial ensemble.
%   V         : d-by-d orthonormal eigenvector matrix.
%   lambda    : d-by-1 modal drift coefficients.
%   mu        : d-by-m modal diffusion coefficients.
%   dt        : time-step size.
%   dW        : m-by-N-by-Nt Brownian increments (already scaled by sqrt(dt)).
%
%   Output
%   ------
%   X : d-by-N-by-(Nt+1) exact ensemble.

    [d,N] = size(X0);
    Nt = size(dW,3);

    Z = V'*X0;
    drift = lambda-0.5*sum(mu.^2,2);

    X = zeros(d,N,Nt+1);
    X(:,:,1) = X0;
    for n = 1:Nt
        Z = Z.*exp(drift*dt+mu*dW(:,:,n));
        X(:,:,n+1) = V*Z;
    end

end
