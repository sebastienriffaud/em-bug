function X = bug_standard(problem,dt,N,dW,r,seed)
%BUG_STANDARD Integrate the ensemble with the standard BUG integrator.
%   X = BUG_STANDARD(problem,dt,N,dW,r,seed) returns the rank-r ensemble
%   trajectory produced by the BUG integrator.
%
%   Inputs
%   ------
%   problem : "ou" or "multiplicative".
%   dt      : positive time-step size.
%   N       : number of particles.
%   dW      : m-by-N-by-Nt Brownian increments (already scaled by sqrt(dt)).
%   r       : prescribed rank.
%   seed    : reserved input for deterministic initialization.
%
%   Output
%   ------
%   X : d-by-N-by-(Nt+1), X(:,:,k) approximates the ensemble at time t(k).

    [d,~,t,Nt,X0,F] = setup(problem,dt,N,seed);

    % Initialization.
    [U,Y] = truncate(X0,r);
    [V,~] = qr(Y,0);

    X = zeros(d,N,Nt+1);
    X(:,:,1) = U*Y';
    for n = 1:Nt
        Xold = U*Y';
        [drift,noise] = F(t(n),Xold,[],dW(:,:,n));

        % Physical basis augmentation.
        Uhat = orth([U (dt*drift+noise)*V]);

        % Stochastic update.
        Yhat = (Xold+dt*drift+noise)'*Uhat;

        % Truncation step.
        [Phi,Sigma,Psi] = svd(Yhat,'econ');
        U = Uhat*Psi(:,1:r);
        Y = Phi(:,1:r)*Sigma(1:r,1:r);
        V = Phi(:,1:r);

        X(:,:,n+1) = U*Y';
    end

end
