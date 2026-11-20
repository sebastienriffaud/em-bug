function X = bug_closed(problem,dt,N,dW,r,seed)
%BUG_CLOSED Integrate the ensemble with the closed-moment EM--BUG integrator.
%   X = BUG_CLOSED(problem,dt,N,dW,r,seed) returns the rank-r ensemble
%   trajectory produced by the closed-moment EM--BUG integrator.
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

    [d,~,t,Nt,X0,F,~,~,W,A] = setup(problem,dt,N,seed);

    % Initialization.
    [U,Y] = truncate(X0,r);
    D = (Y'*Y)/N;
    M = eye(d)+dt*A;

    X = zeros(d,N,Nt+1);
    X(:,:,1) = U*Y';
    for n = 1:Nt
        Xold = U*Y';
        Cold = U*D*U';
        [drift,noise,G,H] = F(t(n),Xold,Cold,dW(:,:,n));

        % Physical basis augmentation.
        Uhat = orth([M*U W(U)]);

        % Second-moment update.
        S = Cold+dt*G+dt^2*H;
        Dhat = Uhat'*S*Uhat;
        Dhat = (Dhat+Dhat')/2; % Enforce symmetry before diagonalization (roundoff protection).

        % Truncation step.
        [Phi,lambda] = eig(Dhat,'vector');
        [~,ind] = sort(lambda,'descend');
        U = Uhat*Phi(:,ind(1:r));
        D = diag(lambda(ind(1:r)));

        % Stochastic update.
        Y = (Xold+dt*drift+noise)'*U;
        
        X(:,:,n+1) = U*Y';
    end

end
