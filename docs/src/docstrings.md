# Docstrings
```@meta
CurrentModule = ExpandLGF
```
There are two ways to use this package for LGF evaluation. 
 - The high level interface performs LGF evaluations by automatically switching from quadrature to a series expansion whenever it is permissible within a given error tolerance. The user specifies the tolerance and the number of series terms used. This functionality is exported by the package.
 - The low level interface provides direct access to expansion coefficients, function generation, and error estimates for each type of series. This functionality is not exported.

## High level interface
### Stencil Types
```@docs
Stencil
SplitStencil
FullStencil
MehrstellenStencil
Mehrstellen4
Mehrstellen6
LeftMehrstellen4
LeftMehrstellen6
StandardDifference2D
StandardDifference3D
```
### LGF Evaluation
```@docs
EvaluateLGF
EvaluateLGF(n)
NearFieldEvaluator
NearFieldSplit
NearFieldFull
NearField
FarFieldExpansion
FarField
```

## Low level interface
### Expansion coefficients
```@docs
inner_integral_expansion_coeffs
outer_integral_expansion_coeffs
outer_integral_2d_expansion_coeffs
outer_integral_3d_expansion_coeffs
graded_series_expansion
far_field_expansion_coeffs
martinsson_rodin_a2q
```

### Callable expansions
```@docs
inner_integral_expansion
outer_integral_expansion
outer_integral_2d_expansion
outer_integral_3d_expansion
far_field_expansion
far_field_2d_expansion
far_field_3d_expansion
```

### Error estimates
```@docs
inner_integral_errors
outer_integral_errors
far_field_errors
```