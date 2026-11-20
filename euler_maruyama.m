function X = euler_maruyama(problem,dt,N,dW,seed)
%EULER_MARUYAMA Integrate the ensemble by the Euler--Maruyama method.
%   X = EULER_MARUYAMA(problem,dt,N,dW,seed) returns the ensemble
%   trajectory produced by the Euler--Maruyama method.
%
%   Inputs
%   ------
%   problem : "ou" or "multiplicative".
%   dt      : positive time-step size.
%   N       : number of particles.
%   dW      : m-by-N-by-Nt Brownian increments (already scaled by sqrt(dt)).
%   seed    : reserved input for deterministic initialization.
%
%   Output
%   ------
%   X : d-by-N-by-(Nt+1), X(:,:,k) approximates the ensemble at time t(k).

    [d,~,t,Nt,X0,F] = setup(problem,dt,N,seed);

    X = zeros(d,N,Nt+1);
    X(:,:,1) = X0;
    for k = 1:Nt
        Xold = X(:,:,k);
        [drift,noise] = F(t(k),Xold,[],dW(:,:,k));
        X(:,:,k+1) = Xold+dt*drift+noise;
    end

end
