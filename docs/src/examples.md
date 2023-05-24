### Generate C code for far-field LGF Expansion

```julia
using ExpandLGF
using Symbolics

stencil = Mehrstellen6()
terms = 3
expansion_coefficients = ExpandLGF.far_field_expansion_coeffs(stencil, terms)
code = ExpandLGF.far_field_3d_expansion(expansion_coefficients, Symbolics.CTarget())
```
The block above generates C code that evaluates a far-field expansion for the sixth-order Mehrstellen stencil. The length of the expansion is determined by `terms `, which indicates the number of terms $m$ in the summation
```math
G(n) = \frac{1}{4\pi |\bm{n}|} + \frac{1}{\pi}\sum_{k=1}^m \frac{g_k(\bm{n})}{|\bm{n}|^{6k+1}} + \mathcal{O}\left( |\bm{n}|^{-(2m + 3)} \right).
```
The expansion computed above will contain the analytical Green's function and the corrections $g_1(\bm{n})$, $g_2(\bm{n})$, and $g_3(\bm{n})$, which provides an accuracy of $\mathcal{O}(n^{-9})$. The output `code` is the string
```cpp
#include <math.h>
void lgf(double* G, const double n1, const double n2, const double n3, const double n) {
  G[0] = (-0.018229166666666668 * pow(n1, 12) + -0.018229166666666668 * pow(n2, 12) + -0.018229166666666668 * pow(n3, 12) + -0.5625 * (n1 * n1) * pow(n2, 10) + 8.088541666666666 * pow(n1, 6) * pow(n2, 6) + 3.5 * pow(n1, 4) * pow(n2, 8) + 8.088541666666666 * pow(n1, 6) * pow(n3, 6) + -0.5625 * (n1 * n1) * pow(n3, 10) + 3.5 * pow(n1, 8) * pow(n2, 4) + 3.5 * pow(n1, 8) * pow(n3, 4) + -0.5625 * pow(n1, 10) * (n2 * n2) + 3.5 * pow(n1, 4) * pow(n3, 8) + -0.5625 * pow(n1, 10) * (n3 * n3) + 3.5 * pow(n2, 4) * pow(n3, 8) + 3.5 * pow(n2, 8) * pow(n3, 4) + 8.088541666666666 * pow(n2, 6) * pow(n3, 6) + -0.5625 * pow(n2, 10) * (n3 * n3) + -0.5625 * (n2 * n2) * pow(n3, 10) + -15.875 * pow(n1, 6) * pow(n2, 4) * (n3 * n3) + -9.953125 * pow(n1, 8) * (n2 * n2) * (n3 * n3) + -15.875 * pow(n1, 4) * (n2 * n2) * pow(n3, 6) + -15.875 * pow(n1, 6) * (n2 * n2) * pow(n3, 4) + -33.1640625 * pow(n1, 4) * pow(n2, 4) * pow(n3, 4) + -15.875 * pow(n1, 4) * pow(n2, 6) * (n3 * n3) + -9.953125 * (n1 * n1) * (n2 * n2) * pow(n3, 8) + -15.875 * (n1 * n1) * pow(n2, 6) * pow(n3, 4) + -15.875 * (n1 * n1) * pow(n2, 4) * pow(n3, 6) + -9.953125 * (n1 * n1) * pow(n2, 8) * (n3 * n3)) / (M_PI * pow(n, 19)) + 1 / (12.566370614359172 * n * 1);
}
```
Note that because the stencil is sixth order, $g_1(\bm{n}) = g_2(\bm{n}) = 0$. For lower order stencils theses terms are nonzero and will appear in the generated C code.

### Precompute a block of 2D LGF values
```julia
using ExpandLGF
using OffsetArrays

N = 8            # Evaluate entries on [0, N - 1]^2
near_terms = 10  # number of terms in near-field expansions
far_terms = 8    # number of terms used in far-field expansion
rtol = 0         # no relative tolerance on evaluations
atol = 1e-15     # absolute tolerance on evaluations

stencil = StandardDifference2D(4)
ELGF = EvaluateLGF(stencil, near_terms, far_terms, rtol, atol)
    
lgf = OffsetArray(zeros(Float64, (N, N)), 0:N-1, 0:N-1)
for index in ExpandLGF.sorted_indices_2D(N - 1)
    lgf[index...] = ELGF(index)
end

ExpandLGF.symmetrize_block!(lgf)
```
Output:
```
8×8 OffsetArray(::Matrix{Float64}, 0:7, 0:7) with eltype Float64 with indices 0:7×0:7:
  0.0       -0.207425  -0.318808  …  -0.495317  -0.519862
 -0.207425  -0.269189  -0.338713     -0.497503  -0.521472
 -0.318808  -0.338713  -0.37609      -0.503719  -0.526114
 -0.384607  -0.393345  -0.414417     -0.513102  -0.533294
 -0.430685  -0.435572  -0.448589     -0.524612  -0.542364
 -0.466274  -0.469411  -0.478132  …  -0.537315  -0.552685
 -0.495317  -0.497503  -0.503719     -0.550506  -0.563713
 -0.519862  -0.521472  -0.526114     -0.563713  -0.575037
```

The code above evaluates a block of LGF entries for the 2D dimension-split stencil of order four. To take advantage of symmetry, the `sorted_indices_2D` function limits the loop to indices $(i, j)$ with $i \le j$.
Once the computation is complete, the rest of the block is filled using the symmetry condition $G(i, j) = G(j, i)$ in the call to `ExpandLGF.symmetrize_block!`. 

For a 3D stencil, `2D` can be replaced with `3D` in the code and the array `lgf` can be allocated with
```julia
lgf = OffsetArray(zeros(Float64, (N, N, N)), 0:N-1, 0:N-1, 0:N-1)
```
In three dimensions, the symmetry condition reduces the number of evaluations from $N^3$ to roughly $N^3 / 6$.