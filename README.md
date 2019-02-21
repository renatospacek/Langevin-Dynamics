# Equilibrium-LD

A numerical framework written in [Julia] that implements Langevin dynamics, used to thermodynamically sample a system of interacting particles consistent with the canonical ensemble. Implemented features include:

* Pairwise interactions: Lennard-Jones potential.
* Cell list: speeds up the force computation, making it O(N) instead of O(N^2).
* Stochastic integration: Euler-Maruyama.
* Periodic boundary conditions.
* Computation of kinetic energy, potential energy, instantaneous temperature.
* Histogram plot of velocity distribution.

Features to be implemented/known limitations include:
* Generalize it to noncubic domains.
* Implement other integrators/interaction models.
* Computation of more observables.
* Efficiency - the current version is far from optimal.

# Usage
All simulation parameters are defined in 'input.jl'. To run, simply execute 'main.jl'.

[Julia]: http://julialang.org
