## Numerical resolution of wave equation in the temporal regime

This mini-project is a standard application of Finite Element Method for elliptic and hyperbolic PDEs with Dirichlet boundary condition. The first part concerns the resolution of Laplace equation with a non-homogeneous elasticity in 
a rectangular domain while in a second time, one deals with the propagation of a wave in this same inhomogeneous medium. The domain of interest throughout the whole mini-project is a rectangle.

### Stationary resolution
The resolution of the Laplace equation via the Finite Element Method is standard here. After an adequate discretization of the Sobolev space, it amounts to solving a linear equation. The associated matrix is known to be sparse; MATLAB takes advantage of such a structure to accelerate the resolution of the linear equation. 

In a homogeneous medium:

<img src="img/Solexacte.jpg" alt="drawing" width="400"/>

In a nonhomogeneous medium:

<img src="img/sigma_2.jpg" alt="drawing" width="400"/> <img src="img/Stat2.jpg" alt="drawing" width="350"/>

On the left is the elasticity coefficient. On the right is the numerical solution of Laplace equation on the rectangular domain.


### Temporal resolution
When dealing with hyperbolic PDEs, it is common to treat the space and time variable differently. The temporal derivative is discretized via a finite-difference scheme, while the spatial part is discretized thanks to the finite element method. More specifically, the temporal discretization is done with the classical leap-frog scheme for two times derivatives. Concerning the computation of the mass matrix, three different techniques are employed: 
- Standard Quadrature,
- Mass condensation,
- and Cholesky decomposition.

In the illustrations below, we have the case of a domain separated into two homogeneous part; the left-half and right-half of the rectangular domain.

![](img/Film_Exo_3_Q5.gif)

![](img/Film_Exo_3_Q6.gif)


For more details on the results, please refer to the report.
