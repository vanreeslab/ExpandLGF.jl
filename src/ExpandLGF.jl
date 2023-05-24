module ExpandLGF

using QuadGK
using AbstractAlgebra
using Symbolics
using OffsetArrays
using HCubature

include("types.jl")
include("utils.jl")
include("common_stencils.jl")
include("near_field_quadrature.jl")
include("near_field_coefficients.jl")
include("far_field_coefficients.jl")
include("callable_expansions.jl")
include("expansion_errors.jl")
include("high_level_interface.jl")

export StandardDifference2D, StandardDifference3D
export Mehrstellen4, Mehrstellen6, LeftMehrstellen4, LeftMehrstellen6
export NearField, FarField, EvaluateLGF

end # module ExpandLGF

