%% Strong-error convergence with respect to the rank

clear; clc; close all;

problem = "ou";
% problem = "multiplicative";

[dt,N,r] = config(problem);
seed = 1;

[~,m,~,Nt,~,~,exact] = setup(problem,dt,N,seed);
[methods,labels,~,markers] = method_list(problem);

% Generate a shared Brownian path and the exact reference.
rng(seed,'twister');
dW = sqrt(dt)*randn(m,N,Nt);
dZ = randn(m,N,Nt);

Xexact = exact(dW,dZ);

Xem = euler_maruyama(problem,dt,N,dW,seed);
error_em = strong_error(Xexact,Xem);

% Evaluate the two BUG integrators for each rank.
errors = zeros(numel(methods),numel(r));
for i = 1:numel(methods)
    for j = 1:numel(r)
        X = methods{i}(problem,dt,N,dW,r(j),seed);
        errors(i,j) = strong_error(Xexact,X);
    end
end

% Plot strong errors.
figure;
hold on;
semilogy(r,error_em*ones(size(r)),'k--','DisplayName','Euler--Maruyama');
for i = 1:numel(methods)
    semilogy(r,errors(i,:),'Marker',markers(i),'DisplayName',labels(i));
end
xlabel('Rank $r$','Interpreter','latex');
ylabel('Error $e_{\rm strong}$','Interpreter','latex');
legend('Location','northeast','Interpreter','latex');
grid on;
