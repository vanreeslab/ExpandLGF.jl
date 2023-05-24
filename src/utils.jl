"""
    semifactorial(x) -> BigInt

Return the semifactorial `x!! = x*(x-2)*...` as a `BigInt`
"""
function semifactorial(x)::BigInt
    y = x
    out = BigInt(1)
    while y > 0
        out *= y
        y -= 2
    end
    return out
end

"""
    fac(x)

Convenience wrapper for BigInt factorial
"""
fac(x) = factorial(big(x))

"""
    ∂σ(cffs, N)

Return the nth derivative of σ(k) for a 1D finite difference with stencil `cffs`
"""
function ∂σ(cffs, N) 
    if N == 0 || N % 2 != 0
        return zero(eltype(cffs))
    end
    return -2*(-1)^(N÷2)*sum(cj*big(j)^N for (j, cj) in enumerate(cffs))
end

# AbstractAlgebra.jl is not optimized for rapid evaluation of polynomials with inexact types.
# This package uses Symbolics.jl to generate pure Julia functions that serve this purpose.

"Convert an AbstractAlgebra.jl univariate polynomial to a Symbolics.jl expression"
function poly_to_symbolic(p::PolyElem{T}, sym; CoeffType=Float64) where T
    # because p may have degree zero, we do not reduce over the (possibly empty) range 0:degree(p)
    out = Num(0)
    for deg = 0:AbstractAlgebra.degree(p)
        out += CoeffType(coeff(p, deg)) * sym^deg
    end
    return out
end

"Convert an AbstractAlgebra.jl multivariate polynomial to a Symbolics.jl expression"
function poly_to_symbolic(p::MPolyElem{T}, syms; CoeffType=Float64) where T
    # because p may have degree zero, we do not reduce over the (possibly empty) coefficients(p) and exponent_vectors(p)
    out = Num(0)
    for it = 1:length(p)
        out += CoeffType(coeff(p, it)) * prod(syms.^exponent_vector(p, it))
    end
    return out
end

"Return all indices in [0, N]^3 that are in descending order"
function sorted_indices_3D(N) 
    indices = NTuple{3, Int}[]
    for i1 in 0:N, i2 in 0:i1, i3 in 0:i2
        push!(indices, (i1, i2, i3))
    end
    return indices
end

"Return all indices in [0, N]^2 that are in descending order"
function sorted_indices_2D(N) 
    indices = NTuple{2, Int}[]
    for i1 in 0:N, i2 in 0:i1
        push!(indices, (i1, i2))
    end
    return indices
end

"Set all permutations of a given index to the same value."
function symmetrize_block!(data::AbstractArray)
    for I in CartesianIndices(data)
        J = sort(abs.([Tuple(I)...]), rev=true)
        data[I] = data[J...]
    end
    return data
end
