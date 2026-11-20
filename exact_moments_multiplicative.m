function [mean_exact,C_exact] = exact_moments_multiplicative(mean0,C0,V,lambda,mu,T)
%EXACT_MOMENTS_MULTIPLICATIVE Exact moments for the multiplicative problem.
%   [mean_exact,C_exact] = EXACT_MOMENTS_MULTIPLICATIVE(...) returns the
%   exact mean and second moment at final time T.
%
%   Inputs
%   ------
%   mean0, C0 : initial mean and non-centered second moment.
%   V         : d-by-d orthonormal eigenvector matrix.
%   lambda    : d-by-1 modal drift coefficients.
%   mu        : d-by-m modal diffusion coefficients.
%   T         : final time.
%
%   Outputs
%   -------
%   mean_exact : d-by-1 exact mean.
%   C_exact    : d-by-d exact non-centered second moment.

    mean_modal = V'*mean0;
    C_modal = V'*C0*V;

    mean_exact = V*(exp(lambda*T).*mean_modal);
    C_exact = V*(exp((lambda+lambda'+mu*mu')*T).*C_modal)*V';

end
