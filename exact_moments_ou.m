function [mean_exact,C_exact] = exact_moments_ou(mean0,C0,V,lambda,sigma,m,T)
%EXACT_MOMENTS_OU Exact moments for the OU problem.
%   [mean_exact,C_exact] = EXACT_MOMENTS_OU(mean0,C0,V,lambda,...
%       sigma,m,T) returns the exact mean and second moment at final time T.
%
%   Inputs
%   ------
%   mean0, C0 : initial mean and non-centered second moment.
%   V, lambda : orthonormal eigenvectors and eigenvalues of the drift.
%   sigma     : additive-noise amplitude.
%   m         : number of forced eigenmodes.
%   T         : final time.
%
%   Outputs
%   -------
%   mean_exact : d-by-1 exact mean.
%   C_exact    : d-by-d exact non-centered second moment.

    decay = exp(lambda*T);
    mean_exact = V*(decay.*(V'*mean0));

    C_modal = V'*C0*V;
    C_modal = (decay*decay').*C_modal;

    variance = zeros(size(lambda));
    for i = 1:m
        if abs(lambda(i))==0.0
            variance(i) = sigma^2*T;
        else
            variance(i) = sigma^2*expm1(2*lambda(i)*T)/(2*lambda(i));
        end
    end

    C_exact = V*(C_modal+diag(variance))*V';

end
