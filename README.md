# BUG integrators for stochastic differential equations

MATLAB implementation of the **closed-moment Euler--Maruyama BUG (EM--BUG)**
integrator for stochastic differential equations (SDEs).

This code supports the paper:

Sébastien Riffaud,
*Robust BUG dynamical low-rank integrators based on the Euler–Maruyama method
for stochastic differential equations. Part I: a closed-moment approach*, arXiv
preprint arXiv:2611.25152 manuscript (2026),
[doi:10.48550/arXiv.2611.25152](https://doi.org/10.48550/arXiv.2611.25152).

If you use this code, please cite the article above. The code implements the
closed-moment EM--BUG integrator and numerical experiments from Section 5: 
Ornstein--Uhlenbeck process with additive noise and linear equation with
commuting multiplicative noise.

## Contents

`test_rank.m` evaluates the strong-error convergence with respect to the rank.
`test_timestep.m` evaluates the strong-error convergence with respect to the time step.
`test_particles.m` evaluates the mean and second-moment errors with respect to the number of particles.
`config.m` provides default experiment parameters.
`setup.m`, `method_list.m` construct the benchmark and select the integrators.
`bug_standard.m`, `bug_closed.m` compute the low-rank BUG solution.
`euler_maruyama.m` computes the full Euler--Maruyama solution.
`ou.m`, `multiplicative.m` construct the drift, noise, and closed moment operators.
`initial_ou.m`, `initial_multiplicative.m` provide the initial ensembles.
`noise_directions_ou.m`, `noise_directions_multiplicative.m` construct the noise-generated directions.
`exact_solution_ou.m`, `exact_solution_multiplicative.m` provide the exact solutions.
`exact_moments_ou.m`, `exact_moments_multiplicative.m` provide the exact final mean and non-centered second moment.
`truncate.m` compresses a low-rank factorization to a prescribed rank.
`aggregate_increments.m` couples the fine- and coarse-grid Brownian increments.
`strong_error.m` computes the maximum-in-time root-mean-square strong error.


## Requirements

- MATLAB **R2023b or newer**.
- No additional MATLAB toolboxes are required by the supplied code.

## Quick start

1. Open MATLAB and change the current folder to the repository directory,
   or add that directory to the MATLAB path.
2. Run one of the experiment scripts:

   ```
   test_rank       % Strong error versus approximation rank
   test_timestep   % Strong error versus time step
   test_particles  % Mean and second-moment errors versus number of particles
   ```

## Author

Sébastien Riffaud

## License

MIT License. Copyright (c) 2026 Sébastien Riffaud. See [`LICENSE`](LICENSE).
