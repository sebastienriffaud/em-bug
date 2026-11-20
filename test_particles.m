%% Errors in the mean and non-centered second moment with respect to the number of particles

clear; clc; close all;

problem = "ou";
% problem = "multiplicative";

nRep = 5;

[dt,~,~,~,~,N,r] = config(problem);
N = N(:)';
Nmax = max(N);

[~,m,~,Nt,~,~,~,moments] = setup(problem,dt,Nmax,1);
[mean_reference,moment_reference] = moments();
[methods,labels] = method_list(problem);

nM = numel(methods);
nR = numel(r);
nN = numel(N);

mean_exact_rep = zeros(nRep,nN);
C_exact_rep = zeros(nRep,nN);
mean_em_rep = zeros(nRep,nN);
C_em_rep = zeros(nRep,nN);
mean_bug_rep = zeros(nM,nR,nN,nRep);
C_bug_rep = zeros(nM,nR,nN,nRep);

% Independent trials; use the same noise for all methods at each N.
for rep = 1:nRep
    fprintf('%s: repetition %d/%d\n',problem,rep,nRep);
    seed = rep;

    rng(seed,'twister');
    dWmax = sqrt(dt)*randn(m,Nmax,Nt);
    dZmax = randn(m,Nmax,Nt);

    for k = 1:nN
        Nk = N(k);
        dW = dWmax(:,1:Nk,:);
        dZ = dZmax(:,1:Nk,:);

        [~,~,~,~,~,~,exact] = setup(problem,dt,Nk,seed);
        X = exact(dW,dZ);
        Xfinal = X(:,:,end);
        mean_exact_rep(rep,k) = norm(mean(Xfinal,2)-mean_reference);
        C_exact_rep(rep,k) = norm(Xfinal*Xfinal'/Nk-moment_reference,'fro');
        clear X Xfinal;

        X = euler_maruyama(problem,dt,Nk,dW,seed);
        Xfinal = X(:,:,end);
        mean_em_rep(rep,k) = norm(mean(Xfinal,2)-mean_reference);
        C_em_rep(rep,k) = norm(Xfinal*Xfinal'/Nk-moment_reference,'fro');
        clear X Xfinal;

        % Evaluate the two BUG integrators at each rank and number of particles.
        for i = 1:nM
            for j = 1:nR
                X = methods{i}(problem,dt,Nk,dW,r(j),seed);
                Xfinal = X(:,:,end);
                mean_bug_rep(i,j,k,rep) = norm(mean(Xfinal,2)-mean_reference);
                C_bug_rep(i,j,k,rep) = norm(Xfinal*Xfinal'/Nk-moment_reference,'fro');
                clear X Xfinal;
            end
        end
    end
end

% Root-mean-square error across the repetitions.
mean_exact_mc = sqrt(mean(mean_exact_rep.^2,1));
C_exact_mc = sqrt(mean(C_exact_rep.^2,1));
mean_em = sqrt(mean(mean_em_rep.^2,1));
C_em = sqrt(mean(C_em_rep.^2,1));
mean_errors = sqrt(mean(mean_bug_rep.^2,4));
C_errors = sqrt(mean(C_bug_rep.^2,4));

reference_mean = sqrt(mean(N.*mean_exact_mc.^2))./sqrt(N);
reference_C = sqrt(mean(N.*C_exact_mc.^2))./sqrt(N);

% Plot errors in the mean and non-centered second moment.
for i = 1:nM
    figure('Name',labels(i),'Position',[600 100 1120 420]);

    subplot(1,2,1);
    for j = 1:nR
        loglog(N,squeeze(mean_errors(i,j,:)),'o-','DisplayName',sprintf('Rank %d',r(j)));
        hold on;
    end
    loglog(N,mean_em,'k--','DisplayName','Euler--Maruyama');
    loglog(N,mean_exact_mc,'k--','Marker','s','DisplayName','Exact Monte Carlo');
    loglog(N,reference_mean,'k:','DisplayName','$O(N^{-1/2})$');
    xlabel('Number of particles $N$','Interpreter','latex');
    ylabel('Error $e_{\rm mean}$','Interpreter','latex');
    title('Mean','Interpreter','latex');
    legend('Location','southwest','Interpreter','latex');
    grid on;

    subplot(1,2,2);
    for j = 1:nR
        loglog(N,squeeze(C_errors(i,j,:)),'o-','DisplayName',sprintf('Rank %d',r(j)));
        hold on;
    end
    loglog(N,C_em,'k--','DisplayName','Euler--Maruyama');
    loglog(N,C_exact_mc,'k--','Marker','s','DisplayName','Exact Monte Carlo');
    loglog(N,reference_C,'k:','DisplayName','$O(N^{-1/2})$');
    xlabel('Number of particles $N$','Interpreter','latex');
    ylabel('Error $e_{\rm moment}$','Interpreter','latex');
    title('Non-centered second moment','Interpreter','latex');
    legend('Location','southwest','Interpreter','latex');
    grid on;
end
