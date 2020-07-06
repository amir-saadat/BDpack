Brownain Dynamics (BD) Simulation
=================================

Brownian dynamics (BD) is an accurate and computationally efficient mesoscale simulation technique used to study the dynamics and material properties of synthetic and biological polymeric solutions with different architectures, under both equilibrium and nonequilibrium (with the presence of an external field) conditions. A common model used to describe the physics of a polymer is the bead-spring model, where beads resemble the centers of hydrodynamic resistance and they are connected by a network of springs, where the spring tension is generally related nonlinearly to the spring extension. In comparison to a fully resolved molecular dynamics model, this bead-spring model abstracts away the fine (unnecessary) details of real macromoleculaes. 

The stochastic differential equation governing the positions of ![equation](http://latex.codecogs.com/gif.latex?\inline&space;N_\mathrm{b}) beads in a bead-spring chain is [[&Ouml;ttinger (1989)]](http://www.springer.com/us/book/9783540583530),

![equation](http://latex.codecogs.com/gif.latex?\dpi{100}&space;\small&space;\text{d}\boldsymbol{r}_{\nu}=\left[Pe&space;\boldsymbol{\kappa}&space;\cdot\boldsymbol{r}_{\nu}&space;&plus;&space;\frac{1}{4}\sum_{\mu=1}^{N_\mathrm{b}}\nabla_{\mu}\cdot&space;\mathbf{D_{\mu\nu}}&space;&plus;\frac{1}{4}&space;\sum_{\mu=1}^{N_\mathrm{b}}\mathbf{D}_{\nu\mu}\cdot&space;\boldsymbol{F}_{\mu}\right]\mathrm{d}t&space;&plus;&space;\frac{1}{\sqrt{2}}&space;\sum_{\mu=1}^{N_\mathrm{b}}\mathbf{C}_{\nu\mu}\cdot&space;\mathrm{d}\boldsymbol{W}_{\mu},)

where d*t* is the time step, *Pe* is the Peclet number, representing the strength of the flow relative to diffusive motion, ![equation](http://latex.codecogs.com/gif.latex?\inline&space;\dpi{100}&space;\small&space;\boldsymbol{\mathrm{\kappa}}) is the transposed, gradient of the imposed velocity field, and ![equation](http://latex.codecogs.com/gif.latex?\inline&space;\dpi{100}&space;\small&space;\mathrm{d}\boldsymbol{W}) is an increment of a 3![equation](http://latex.codecogs.com/gif.latex?\inline&space;N_\mathrm{b}) dimensional Wiener process. The net non-hydrodynamic, non-Brownian force acting on each bead is

![equation](http://latex.codecogs.com/gif.latex?\dpi{100}&space;\small&space;\boldsymbol{F}&space;=&space;\boldsymbol{F}^{\phi}&space;&plus;\boldsymbol{F}^\mathrm{ext})

where ![equation](http://latex.codecogs.com/gif.latex?\inline&space;\dpi{100}&space;\small&space;\boldsymbol{F}^{\phi}) is all conservative forces, i.e., the net entropic spring force, the excluded volume force preventing beads from intersecting, and bending forces. ![equation](http://latex.codecogs.com/gif.latex?\inline&space;\boldsymbol{F}^{\text{ext}}) is any external force that may be applied. The hydrodynamic interactions (HI) between beads are captured by  ![equation](http://latex.codecogs.com/gif.latex?\inline&space;\mathbf{D}_{\nu\mu}) which relates forces acting on a bead at position ![equation](http://latex.codecogs.com/gif.latex?\inline&space;\boldsymbol{x}_\mu) to a velocity perturbation on a bead at location ![equation](http://latex.codecogs.com/gif.latex?\inline&space;\boldsymbol{x}_\nu), i.e., 

![equation](http://latex.codecogs.com/gif.latex?\boldsymbol{v}^{\prime}(\boldsymbol{x}_\nu)&space;=&space;\mathbf{D}_{\nu\mu}\cdot&space;\boldsymbol{F}(\boldsymbol{x}_\mu).)

The form of ![equation](http://latex.codecogs.com/gif.latex?\inline&space;\mathbf{D}) is available analytically for polymers in simple domains such as a polymer in free space or a polymer near a plane wall, and for more complex domains this may be determined numerically. As a result of the fluctuation-dissipation theorem, we require that 

![equation](http://latex.codecogs.com/gif.latex?\mathbf{D}&space;=&space;\mathbf{C}\cdot&space;\mathbf{C}^\mathrm{T}.)

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
