# LGFOneUnbounded
A small header-only C++ library for calculating Lattice Green's Functions (LGFs) on a Cartesian grid with one unbounded dimension and two periodic dimensions. 

The header `lgf_one_unbounded.hpp` includes the Green's function and Fourier symbol for:
 - second, fourth, sixth, and eighth order dimension-split stencils
 - fourth and sixth order Mehrstellen stencils

### Functionality
This library does not directly evaluate the LGF `G(n0, n1, n2)` that is unbounded and symmetric in `n0` and periodic in `n1, n2` with period `N1, N2`. Instead, the code evaluates the partially transformed LGF `G(n, k1, k2)` which has been expanded in a Fourier series along the `n1` and `n2` directions. To reconstruct the full LGF a discrete Fourier transform is necessary, using wavenumbers `k1 = {2*M_PI*n1/N1, 0 <= n < N1}` along dimension one and likewise for `k2` along dimension two.

Some other considerations:
 - Greens functions are for the negative Laplacian operator.
 - All functions are provided in the namespace `LGFOneUnbounded`
 - Calculations are designed for `double` precision, but this can be reduced to `float` precision by redefining the type `LGFOneUnbounded::real_t`.
 - The library assumes that the LGF is defined on a grid with unit grid spacing, so when working on a grid with spacing `h != 1.0` the result should be rescaled.
 - For 2D domains with one unbounded direction, set `k2 = 0`.

### Header-only usage
Copy `lgf_one_unbounded.hpp` from the top level directory to your source tree, and `#include "lgf_one_unbounded.hpp"` in your source code. 

```cpp
#include "lgf_one_unbounded.hpp"

using LGFOneUnbounded;

int n = 10;
double k1 = M_PI, k2 = M_PI;

lgf2(n, k1, k2); // second order dimension-split stencil
lgf4(n, k1, k2); // fourth order dimension-split stencil
lgf6(n, k1, k2); // sixth order dimension-split stencil
lgf8(n, k1, k2); // eight order dimension-split stencil

meh4_left(n, k1, k2); // fourth order Mehrstellen stencil, left hand side only
meh4_full(n, k1, k2); // fourth order Mehrstellen stencil, combined left and right stencils
meh6_left(n, k1, k2); // sixth order Mehrstellen stencil, left hand side only
meh6_full(n, k1, k2); // sixth order Mehrstellen stencil, combined left and right stencils
```

### Testing
To build the tests you will need a C++ compiler, the Googletest library, and GNU make. For a standard Googletest installation, setting `GTEST_DIR` to the installation directory is sufficient to build the tests:
```
GTEST_DIR=/path/to/googletest make
./lgf_test
```
Otherwise the path to the Googletest headers and libraries can be specificed directly, along with the names of the Googletest libraries -- see the Makefile for more detail.
