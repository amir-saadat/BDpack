Brownain Dynamics (BD) Simulation
=================================

Brownian dynamics (BD) is an accurate and computationally efficient mesoscale simulation technique used to study the dynamics and material properties of synthetic and biological polymeric solutions with different architectures, under both equilibrium and nonequilibrium (with the presence of an external field) conditions. A common model used to describe the physics of a polymer is the bead-spring model, where beads resemble the centers of hydrodynamic resistance and they are connected by a network of springs, where the spring tension is generally related nonlinearly to the spring extension. In comparison to a fully resolved molecular dynamics model, this bead-spring model abstracts away the fine (unnecessary) details of real macromoleculaes. 

The stochastic differential equation governing the positions of $N_\mathrm{b}$ beads in a bead-spring chain is [[&Ouml;ttinger (1989)]](http://www.springer.com/us/book/9783540583530),

$$\text{d}\boldsymbol{r}_{\nu}=\left[Pe\boldsymbol{\kappa}\cdot\boldsymbol{r}_{\nu}+\frac{1}{4}\sum_{\mu=1}^{N_\mathrm{b}}\nabla_{\mu}\cdot\mathbf{D_{\mu\nu}}+\frac{1}{4}\sum_{\mu=1}^{N_\mathrm{b}}\mathbf{D}_{\nu\mu}\cdot\boldsymbol{F}_{\mu}\right]\mathrm{d}t+\frac{1}{\sqrt{2}}\sum_{\mu=1}^{N_\mathrm{b}}\mathbf{C}_{\nu\mu}\cdot\mathrm{d}\boldsymbol{W}_{\mu}$$

where dt is the time step, Pe is the Peclet number, representing the strength of the flow relative to diffusive motion, $\boldsymbol{\kappa}$ is the transposed, gradient of the imposed velocity field, and $\boldsymbol{W}$ is an increment of a $3N_\mathrm{b}$ dimensional Wiener process. The net non-hydrodynamic, non-Brownian force acting on each bead is

$$\boldsymbol{F}=\boldsymbol{F}^{\phi}+\boldsymbol{F}^\mathrm{ext}$$

<!-- where $\boldsymbol{F}^{\phi}$ is all conservative forces, i.e., the net entropic spring force, the excluded volume force preventing beads from intersecting, and bending forces. $\boldsymbol{F}^{\text{ext}}$ is any external force that may be applied. The hydrodynamic interactions (HI) between beads are captured by  $\mathbf{D}_{\nu\mu}$ which relates forces acting on a bead at position $\boldsymbol{x}_\mu$ to a velocity perturbation on a bead at location $\boldsymbol{x}_\nu$, i.e.,  -->

$$\boldsymbol{v}^{\prime}(\boldsymbol{x}_\nu)=\mathbf{D}_{\nu\mu}\cdot\boldsymbol{F}(\boldsymbol{x}_\mu).$$

The form of ![equation](http://latex.codecogs.com/gif.latex?\inline&space;\mathbf{D}) is available analytically for polymers in simple domains such as a polymer in free space or a polymer near a plane wall, and for more complex domains this may be determined numerically. As a result of the fluctuation-dissipation theorem, we require that 

$$\mathbf{C}\cdot\mathbf{C}^\mathrm{T}.$$

The diffusivity matrix ![equation](http://latex.codecogs.com/gif.latex?\inline&space;\mathbf{D}) is a function of the current bead coordinates ![equation](http://latex.codecogs.com/gif.latex?\inline&space;\boldsymbol{r}), and thus must be recalculated at every time step of the simulation. As a result, the matrix ![equation](http://latex.codecogs.com/gif.latex?\inline&space;\mathbf{C}) must also be recalculated through a decomposition of ![equation](http://latex.codecogs.com/gif.latex?\inline&space;\mathbf{D}) as given by the equation given above.

* [Documentation](https://github.com/amir-saadat/BDpack/wiki/Documentation)
  + [Installation and Running](https://github.com/amir-saadat/BDpack/wiki/Installation-and-Running)
  + User Inputs
    - [Dilute Linear Polymers](https://github.com/amir-saadat/BDpack/wiki/User-Inputs-(DiluteLinear))
    - [Dilute Linear Polymers](https://github.com/amir-saadat/BDpack/wiki/User-Inputs-(DiluteComb))
    - [Semi-dilute Linear Polymers](https://github.com/amir-saadat/BDpack/wiki/User-Inputs-(SemidiluteLinear))
  + Tutorials
    - [Dilute Linear Polymers](https://github.com/amir-saadat/BDpack/wiki/Tutorials-(DiluteLinear))
    - [Dilute Linear Polymers](https://github.com/amir-saadat/BDpack/wiki/Tutorials-(DiluteComb))
    - [Semi-dilute Linear Polymers](https://github.com/amir-saadat/BDpack/wiki/Tutorials-(SemidiluteLinear))
* [Citation](https://github.com/amir-saadat/BDpack/wiki/Citation)
* [Publications](https://github.com/amir-saadat/BDpack/wiki/Publications)
* [License](https://github.com/amir-saadat/BDpack/wiki/License)
