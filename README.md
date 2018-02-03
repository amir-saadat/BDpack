BDpack
======
BDpack is a package to numerically calculate the configurational evolution of polymeric solution using Brownian dynamics simulation of bead-spring micro-mechanical model. The algorithm used in BDpack incorporates high-fidelity and computationally efficient calculation of hydrodynamic interactions (HI) and excluded volume (EV) forces. BDpack is written in parallel by employing message passing interface (MPI) for distributed-memory-architectures.

The aim of this project is to write the codes self descriptive, documented, yet very efficient to enable its applicability and development. To this end, the codes are written in modular and object-oriented fashion using Fortran.

Brownian dynamics (BD) is an accurate and computationally efficient mesoscale simulation technique used to study the dynamics and material properties of synthetic and biological polymeric solutions with different architectures, under both equilibrium and nonequilibrium (with the presence of an external field) conditions. A common model used to describe the physics of a polymer is the bead-spring model, where beads resemble the centers of hydrodynamic resistance and they are connected by a network of springs, where the spring tension is generally related nonlinearly to the spring extension. In comparison to a fully resolved molecular dynamics model, this bead-spring model abstracts away the fine (unnecessary) details of real macromoleculaes. 

The stochastic differential equation governing the positions of $N_\mathrm{b}$ beads in a bead-spring chain is 

<img src="https://github.com/amir-saadat/BDpack/blob/master/img.png" width="600">

* [Documentation](https://github.com/amir-saadat/BDpack/wiki/Documentation)
  + [Installation and Running](https://github.com/amir-saadat/BDpack/wiki/Installation-and-Running)
  + User Inputs
    - [Dilute](https://github.com/amir-saadat/BDpack/wiki/User-Inputs-(Dilute))
    - [Semi-dilute](https://github.com/amir-saadat/BDpack/wiki/User-Inputs-(Semidilute))
  + Tutorials
    - [Dilute](https://github.com/amir-saadat/BDpack/wiki/Tutorials-(Dilute))
    - [Semi-dilute](https://github.com/amir-saadat/BDpack/wiki/Tutorials-(Semidilute))
* [Citation](https://github.com/amir-saadat/BDpack/wiki/Citation)
* [Publications](https://github.com/amir-saadat/BDpack/wiki/Publications)
* [License](https://github.com/amir-saadat/BDpack/wiki/License)
