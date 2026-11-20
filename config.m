function [dt,N,rank_r,timestep_r,timestep_dt,particle_N,particle_r] = config(problem)
%CONFIG Default parameters for the numerical experiments.
%   [dt,N,rank_r,timestep_r,timestep_dt,particle_N,particle_r] = ...
%       CONFIG(problem) returns experiment settings for the selected SDE.
%
%   Input
%   -----
%   problem : "ou" or "multiplicative".
%
%   Outputs
%   -------
%   dt          : time step for rank and particle experiments.
%   N           : number of particles for rank and time-step experiments.
%   rank_r      : ranks used in TEST_RANK.
%   timestep_r  : ranks used in TEST_TIMESTEP.
%   timestep_dt : time steps used in TEST_TIMESTEP.
%   particle_N  : number of particles used in TEST_PARTICLES.
%   particle_r  : ranks used in TEST_PARTICLES.

    switch problem
        % OU problem parameters.
        case "ou"
            dt = 1/128;
            N = 5000;
            rank_r = 1:20;
            timestep_r = [13 14 15];
            timestep_dt = 1./[1 2 4 8 16 32 64 128];
            particle_N = [500 1000 2000 5000];
            particle_r = [5 10 15];

        % Multiplicative-noise problem parameters.
        case "multiplicative"
            dt = 1/128;
            N = 5000;
            rank_r = 1:20;
            timestep_r = [5 7 10 15];
            timestep_dt = 1./[1 2 4 8 16 32 64 128];
            particle_N = [500 1000 2000 5000];
            particle_r = [5 7 15];

        otherwise
            error('config:UnknownTestcase','Unknown testcase "%s".',problem)
    end

end
