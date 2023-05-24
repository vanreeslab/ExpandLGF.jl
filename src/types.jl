"An abstract type for all finite difference stencils"
abstract type Stencil{N, T} end

"A dimension-split finite difference stencil"
struct SplitStencil{N, T} <: Stencil{N, T}
    coefficients::Vector{T}
end
width(ss::SplitStencil) = length(ss.coefficients)
coefficients(ss::SplitStencil) = ss.coefficients

function Base.show(io::IO, stencil::SplitStencil{N, T}) where {N, T}
    print(io, "$(N)D SplitStencil with coefficients $(stencil.coefficients)")
end

"A general symmetric finite difference stencil in 2D or 3D"
struct FullStencil{N, T} <: Stencil{N, T}
    width::Int
    indices::Vector{NTuple{N, Int}}
    values::Vector{T}
end
width(fs::FullStencil) = fs.width
indices(fs::FullStencil) = fs.indices
values(fs::FullStencil) = fs.values

"""
    FullStencil(coefficients)

Construct a `FullStencil` by providing a dictionary of coefficiencts. 

The dictionary should take advantage of symmetry by providing each unique coefficient
only once. The corresponding index should be sorted from low to high. As an example,
```
Dict(
    (0,0,0) => -64 // 15,
    (0,0,1) => +7 // 15,
    (0,1,1) => +1 // 10,
    (1,1,1) => +1 // 30,
)
```
represents a 27-point stencil with given coefficients for the center, face, edge, and corner points.
"""
function FullStencil(coefficients::Dict{NTuple{N, Int}, T}) where {N, T}
    @assert(N in (2, 3), "Only 2D and 3D stencils are supported")
    assert_dict_format_(coefficients)
    stencil_width = stencil_width_from_dict_(coefficients)
    indices, values = indices_and_values_from_dict_(coefficients, stencil_width)

    # final sanity checks on the input stencil
    @assert(sum(values) == zero(T), "Input stencil is not consistent.")

    return FullStencil{N, T}(stencil_width, indices, values)
end

"A Mehrstellen stencil with both left and right finite difference stencils"
struct MehrstellenStencil{N, T} <: ExpandLGF.Stencil{N, T}
    left_stencil::FullStencil{N, T}
    right_stencil::FullStencil{N, T}
end
left_stencil(ms::MehrstellenStencil) = ms.left_stencil
right_stencil(ms::MehrstellenStencil) = ms.right_stencil

"""
    FullStencil(coefficients)

Construct a `FullStencil` by providing a dictionary of coefficiencts for the left and right operators.

See `FullStencil` for the expected format of the `coefficients` argument.
"""
function MehrstellenStencil(lhs_coefficients::Dict{NTuple{N, Int}, T}, rhs_coefficients::Dict{NTuple{N, Int}, T}) where {N, T}
    # RHS stencil does not need to be consistent, so we do not directly use the FullStencil constructor
    assert_dict_format_(rhs_coefficients)
    rhs_stencil_width = stencil_width_from_dict_(rhs_coefficients)
    rhs_indices, rhs_values = indices_and_values_from_dict_(rhs_coefficients, rhs_stencil_width)
    rhs_stencil = FullStencil{N, T}(rhs_stencil_width, rhs_indices, rhs_values)
    lhs_stencil = FullStencil(lhs_coefficients)
    return MehrstellenStencil{N, T}(lhs_stencil, rhs_stencil) 
end


function assert_dict_format_(coefficients::Dict{NTuple{N, Int}, T}) where {N, T}
    indices_are_positive = all(all(index .>= zero(T)) for index in keys(coefficients))
    indices_are_ordered = all(sort([index...]) == [index...] for index in keys(coefficients))
    @assert(indices_are_positive && indices_are_ordered)
end

function stencil_width_from_dict_(coefficients)
    stencil_width = 0
    for cartesian_index in keys(coefficients)   
        for index in cartesian_index
            stencil_width = max(stencil_width, abs(index))
        end
    end
    return stencil_width
end

function indices_and_values_from_dict_(coefficients::Dict{NTuple{N, Int}, T}, stencil_width::Int) where {N, T}
    # expand the symmetrized input into a full list of index / value pairs
    indices = NTuple{N, Int}[]
    values = T[]
    stencil_range = -stencil_width:stencil_width
    product_range = Iterators.product(fill(stencil_range, N)...)
    for I in product_range
        symmetrized_index = Tuple(sort(abs.([I...])))
        if symmetrized_index in keys(coefficients)
            push!(indices, I)
            push!(values, coefficients[symmetrized_index...])
        end
    end
    return indices, values
end

### 
# Typing for Internal Computations
###
# AbstractAlgebra.Floats{T} if T <: Float
# QQ if stencil coefficients are <: Rational or Int

DefaultRing(::Type{<:Integer}) = QQ
DefaultRing(::Type{<:Rational}) = QQ
DefaultRing(::Type{<:AbstractFloat}) = RDF

function DefaultUnivariatePolyRing(::Stencil{N, T}, s::String = "x") where {N, T}
    return PolynomialRing(DefaultRing(T), s)
end

function DefaultMultivariatePolyRing(::Stencil{N, T}, s::String = "x") where {N, T}
    return PolynomialRing(DefaultRing(T), [s*"$(i)" for i in 1:N])
end

function DefaultMultivariateSeriesRing(::Stencil{N, T}, precision::Int, s::String = "x") where {N, T}
    return PowerSeriesRing(DefaultRing(T), fill(precision, N), [s*"$(i)" for i in 1:N])
end