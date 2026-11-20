function [drift,noise,G,H] = multiplicative(~,X,C,V,lambda,mu,dW)
%MULTIPLICATIVE Multiplicative problem operators.
%   [drift,noise,G,H] = MULTIPLICATIVE(t,X,C,V,lambda,mu,dW) evaluates the
%   drift, noise, and moment operators.
%
%   Inputs
%   ------
%   t          : time (unused; this model is autonomous).
%   X          : d-by-N ensemble.
%   C          : d-by-d non-centered second moment.
%   V          : d-by-d modal basis.
%   lambda, mu : modal drift and diffusion coefficients.
%   dW         : m-by-N Brownian increments.
%
%   Outputs
%   -------
%   drift, noise : d-by-N Euler--Maruyama updates.
%   G, H         : d-by-d moment-update matrices.

    Z = V'*X;
    drift = V*(lambda.*Z);
    noise = V*(Z.*(mu*dW));

    if nargout > 2
        Cmodal = V'*C*V;
        Gmodal = (lambda+lambda').*Cmodal+(mu*mu').*Cmodal;
        Hmodal = (lambda*lambda').*Cmodal;
        G = V*Gmodal*V';
        H = V*Hmodal*V';
    end

end
