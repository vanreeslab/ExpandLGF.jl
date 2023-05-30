# ExpandLGF.jl

ExpandLGF.jl provides fast and accurate evaluations of the Lattice Green's Functions (LGFs) for high-order finite difference discretizations of Poisson's equation.

## Background
Lattice Green's Functions (LGFs) are fundamental solutions to discretized linear operators that are frequently used in the numerical solution to elliptic PDEs. As such, they are a useful tool for solving Poisson's equation on domains with one or more free-space boundary conditions, which appear frequently in simulations of incompressible fluid flows (see [ViscousFlow.jl](https://github.com/JuliaIBPM/ViscousFlow.jl) for one such example). 

Generally LGFs are defined through a Fourier integral with a singular and oscillatory integrand. For the standard second-order centered discretization of the Laplacian, analytical techniqes can be used to greatly simplify the 2D or 3D integral into a 1D form which is amenable to numerical quadrature. This package uses a combination of symbolic computation and numerical quadrature to extend these simplifications to higher-order discretizations. 

For interactions between distant sources, a far-field asymptotic expansion of the LGF can bypass numerical integration altogether. This package implements an existing algorithm from [Martinsson and Rodin](https://doi.org/10.1098/rspa.2002.0985) to calculate these far-field expansions. The results are available either as compiled Julia functions or C code.

Finally, the LGF concept can be extended to domains with a combination of periodic and free-space boundary conditions. For domains with one unbounded dimension and one or more periodic dimensions, a C++ header that evaluates the LGFs of several common discretizations is included in the `ccode` folder.