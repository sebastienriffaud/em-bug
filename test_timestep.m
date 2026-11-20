%% Strong-error convergence with respect to the time-step size

clear; clc; close all;

% problem = "ou";
problem = "multiplicative";

[~,N,~,r,dt] = config(problem);
dt_fine = min(dt);
seed = 1;

[~,m,~,Nt_fine,~,~,exact_fine] = setup(problem,dt_fine,N,seed);
[methods,labels] = method_list(problem);

% Generate all Brownian increments on the finest time grid.
rng(seed,'twister');
dWfine = sqrt(dt_fine)*randn(m,N,Nt_fine);
dZfine = randn(m,N,Nt_fine);
Xexact_fine = exact_fine(dWfine,dZfine);

error_em = zeros(size(dt));
errors = zeros(numel(methods),numel(r),numel(dt));

% Aggregate increments consistently for every coarser grid.
for k = 1:numel(dt)
    factor = round(dt(k)/dt_fine);
    dW = aggregate_increments(dWfine,factor);
    Xexact = Xexact_fine(:,:,1:factor:end);

    Xem = euler_maruyama(problem,dt(k),N,dW,seed);
    error_em(k) = strong_error(Xexact,Xem);
    for i = 1:numel(methods)
        for j = 1:numel(r)
            X = methods{i}(problem,dt(k),N,dW,r(j),seed);
            errors(i,j,k) = strong_error(Xexact,X);
        end
    end
end

% Plot strong errors.
reference_one = error_em(end)*dt/dt(end);
reference_half = error_em(end)*sqrt(dt/dt(end));
for i = 1:numel(methods)
    figure('Name',labels(i));
    hold on;
    for j = 1:numel(r)
        loglog(dt,squeeze(errors(i,j,:)),'o-','DisplayName',sprintf('Rank %d',r(j)));
    end
    loglog(dt,error_em,'k--','DisplayName','Euler--Maruyama');
    if problem == "ou"
        loglog(dt,reference_one,'k:','DisplayName','$O(\Delta t)$');
    else
        loglog(dt,reference_half,'k:','DisplayName','$O(\sqrt{\Delta t})$');
    end
    xlabel('Time-step size $\Delta t$','Interpreter','latex','FontSize',13);
    ylabel('Error $e_{\rm strong}$','Interpreter','latex','FontSize',14);
    legend('Location','southeast','Interpreter','latex','FontSize',12);
    xscale("log");
    yscale("log");
    fontsize(scale=1.2);
    grid on;
end
