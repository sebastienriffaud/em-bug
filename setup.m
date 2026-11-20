function [d,m,t,Nt,X0,F,exact,moments,W,A] = setup(problem,dt,N,seed)
%SETUP Construct an SDE benchmark and its reference solutions.
%   [d,m,t,Nt,X0,F,exact,moments,W,A] = SETUP(problem,dt,N,seed) returns
%   the model parameters, the initial ensemble, and function handles.
%
%   Inputs
%   ------
%   problem : "ou" or "multiplicative".
%   dt      : positive time-step size.
%   N       : positive integer, number of particles.
%   seed    : nonnegative integer (reserved input for deterministic initialization).
%
%   Outputs
%   -------
%   d, m    : state dimension and dimension of the noise basis.
%   t, Nt   : time grid and number of time steps.
%   X0      : d-by-N initial ensemble.
%   F       : handle for [drift,noise,G,H] = F(t,X,C,dW).
%   exact   : handle returning the exact ensemble.
%   moments : handle returning the exact final mean and second moment.
%   W       : handle returning the noise-generated basis directions.
%   A       : d-by-d deterministic drift matrix.

    arguments
        problem (1,1) string
        dt (1,1) double {mustBePositive}
        N (1,1) double {mustBeInteger,mustBePositive}
        seed (1,1) double {mustBeInteger,mustBeNonnegative}
    end

    switch problem
        % Additive-noise Ornstein--Uhlenbeck problem.
        case "ou"
            d = 50;
            m = 14;
            T = 1;

            alpha = 1;
            beta = 0.1;
            sigma = 0.08;

            e = ones(d,1);
            L = spdiags([e -2*e e],-1:1,d,d);
            A = -alpha*eye(d)+beta*full(L);
            [V,Lambda] = eig(A);
            lambda = diag(Lambda);

            Vk = V(:,1:m);
            B = sigma*Vk;
            W = @(U) noise_directions_ou(U,Vk);
            [X0,mean0,C0] = initial_ou(d,N,seed);

            F = @(t,X,C,dW) ou(t,X,C,A,B,dW);
            exact = @(dW,dZ) exact_solution_ou(X0,V,lambda,sigma,m,dt,dW,dZ);
            moments = @() exact_moments_ou(mean0,C0,V,lambda,sigma,m,T);

        % Multiplicative-noise problem.
        case "multiplicative"
            d = 50;
            m = 4;
            k = 15;
            T = 1;

            e = ones(d,1);
            L = spdiags([e -2*e e],-1:1,d,d);
            [V,~] = eig(full(L));

            s = linspace(0,1,d)';
            lambda = -0.6-1.2*s;
            mu = zeros(d,m);
            for j = 1:m
                mu(:,j) = 1.2*cos(j*pi*((1:d)'-0.5)/k)/sqrt(m);
            end

            A = V*diag(lambda)*V';
            W = @(U) noise_directions_multiplicative(U,V,mu);
            [X0,mean0,C0] = initial_multiplicative(N,V,k,seed);

            F = @(t,X,C,dW) multiplicative(t,X,C,V,lambda,mu,dW);
            exact = @(dW,dZ) exact_solution_multiplicative(X0,V,lambda,mu,dt,dW);
            moments = @() exact_moments_multiplicative(mean0,C0,V,lambda,mu,T);

        otherwise
            error('setup:UnknownTestcase','Unknown testcase "%s".',problem)
    end

    Nt = round(T/dt);
    t = (0:Nt)*dt;

end
