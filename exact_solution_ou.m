function X = exact_solution_ou(X0,V,lambda,sigma,m,dt,dW,dZ)
%EXACT_SOLUTION_OU Exact ensemble for the OU problem.
%   X = EXACT_SOLUTION_OU(X0,V,lambda,sigma,m,dt,dW,dZ) returns the exact
%   ensemble.
%
%   Inputs
%   ------
%   X0        : d-by-N initial ensemble.
%   V, lambda : eigenvectors and eigenvalues of the drift.
%   sigma     : additive-noise amplitude.
%   m         : number of forced eigenmodes.
%   dt        : time-step size.
%   dW        : m-by-N-by-Nt Brownian increments (already scaled by sqrt(dt)).
%   dZ        : m-by-N-by-Nt independent standard normal variables.
%
%   Output
%   ------
%   X : d-by-N-by-(Nt+1) exact ensemble.

    [d,N] = size(X0);
    Nt = size(dW,3);

    Z = V'*X0;
    decay = exp(lambda*dt);

    coefficient = zeros(m,1);
    bridge = zeros(m,1);
    for i = 1:m
        if abs(lambda(i))==0.0
            variance = dt;
            covariance = dt;
        else
            variance = expm1(2*lambda(i)*dt)/(2*lambda(i));
            covariance = expm1(lambda(i)*dt)/lambda(i);
        end

        coefficient(i) = covariance/dt;
        bridge(i) = sqrt(max(variance-covariance^2/dt,0));
    end

    X = zeros(d,N,Nt+1);
    X(:,:,1) = X0;
    for k = 1:Nt
        Z = decay.*Z;
        Z(1:m,:) = Z(1:m,:)+sigma*(coefficient.*dW(:,:,k)+bridge.*dZ(:,:,k));
        X(:,:,k+1) = V*Z;
    end

end
