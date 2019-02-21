# Langevin-Dynamics

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
All simulation parameters are defined as constants in `input.jl`. To run the code, simply execute `main.jl`. A brief description of the other scrips is below.

* `cell_lists.jl`: Initializes and unwraps neighbor lists to compute the force due to pairwise interactions.
* `hamiltonian.jl`: Computes kinetic/potential energy and the force.
* `initialize.jl`: Generates initial lattice and initial velocities.
* `integrator.jl`: Euler-Maruyama integrator.
* `PBC.jl`: Periodic boundary conditions.
* `histogram.jl`: Normalizes the Gibbs distribution and plots velocity histogram.

# Notes
This code was written as a test code, as this was my first experience with Julia.

I am currently working on a larger project, which implements nonequilibrium Langevin dynamics for a series of background flows, for which this code serves as a basis/prototype. 

All feedback is greatly welcome.

[Julia]: http://julialang.org
