@doc raw"""
    bj = graded_series_expansion(stencil::Stencil{T}, j_max::Int, PR::MPolyRing{T}) where T

Return a power series expansion of the Fourier symbol for a given stencil:
```math
\sigma(\bm{k}) = \sum_{j = 0}^{j_{\mathrm{max}}} b_j(\bm{k}) + \mathcal{O}\left( |\bm{k}|^{2j_{\mathrm{max}} + 2} \right)
```
In the output, `b[j]` contains a homogeneous multivariate polynomial of degree `2j`.

# Arguments
 - `stencil` : the stencil defining the symbol
 - `j_max` : the number of valid entries in the output `bj`
 - `PR` : a multivariate `PolynomialRing` defining the type of the expansion
"""
function graded_series_expansion(stencil::SplitStencil, j_max::Int, PR::MPolyRing) where {N, T}
    b = OffsetArray(zeros(PR, j_max + 1), 0:j_max)
    for j in 1:j_max
        for (i, cff) in enumerate(coefficients(stencil))
            for k in gens(PR)
                b[j] += cff * (2 * (-1)^(j-1) // fac(2*j)) * (i*k)^(2*j)
            end
        end
    end
    return b
end

function graded_series_expansion(stencil::FullStencil{3, T}, j_max::Int, PR::MPolyRing) where {T}
    k = gens(PR)
    b = OffsetArray(zeros(PR, j_max + 1), 0:j_max)
    for β1 = 2*(0:j_max), β2 = 2*(0:j_max), β3 = 2*(0:j_max)
        
        β = (β1, β2, β3)
        j = sum(β) ÷ 2
        if (j > j_max)
            continue
        end

        prefactor = sum(s * prod(α.^β) for (α, s) in zip(indices(stencil), values(stencil)))
        b[j] += (prefactor * (-1)^(j-1) // prod(fac.(β))) * prod(k.^β) 
    end
    return b
end

function graded_series_expansion(stencil::MehrstellenStencil{N, T}, j_max::Int, MPR::MPolyRing) where {N, T}
    # get series expansions for the numerator and denominator
    bl = graded_series_expansion(left_stencil(stencil), j_max, MPR)
    br = graded_series_expansion(right_stencil(stencil), j_max, MPR)

    # convert to multivariate series and divide
    PSR, k_series = DefaultMultivariateSeriesRing(stencil, 2*j_max + 2, "k")
    lhs_series = PSR(+1) * sum(bj(k_series...) for bj in bl)
    rhs_series = PSR(-1) * sum(bj(k_series...) for bj in br) # sign convention
    out_series = divexact(lhs_series, rhs_series)

    # convert back to multivariate polynomials, collected by degree
    k_mpr = gens(MPR)
    b = OffsetArray(zeros(MPR, j_max + 1), 0:j_max)
    for (cff, ev) in zip(AbstractAlgebra.coefficients(out_series), AbstractAlgebra.exponent_vectors(out_series))
        j = sum(ev) ÷ 2
        if (j > j_max)
            continue
        end
        b[j] += cff * prod(k_mpr.^ev)
    end
    return b
end

@doc raw"""
    a2q = martinsson_rodin_a2q(stencil::Stencil, q_max::Int, PR::MPolyRing)

Return the polynomials $a_{2q}(\bm{k})$ with `q <= q_max` for a given stencil.

See Martinnson, Rodin. "Asymptotic Expansions of Lattice Green's Functions". JCP 2022

# Arguments
 - `cffs` : the finite difference coefficients
 - `q_max` : the first `q_max` polynomials are returned
 - `PR` : a vector of 3 generators from a multivariate `PolynomialRing`
"""
function martinsson_rodin_a2q(stencil::Stencil{N, T}, q_max::Int, PR::MPolyRing) where {N, T}

    # Power series expansion of the symbol
    b = graded_series_expansion(stencil, q_max + 1, PR)

    # Martinnson and Rodin's recursion for 1/σ = 1/|x|^2 + sum_q=1^inf a2q / x^(2q+2) 
    PolyType = elem_type(PR)
    c = Dict{Tuple{Int64, Int64}, PolyType}()

    k = gens(PR)
    normk2 = sum(k.^2)

    for j = 2:q_max+1
        c[(0,j)] = -b[j]
    end

    for i1 = 3:q_max+1
        for i2 = 1:i1-2
            j = i1 + i2 
            q = i2
            c[(q,j)] = -b[j-2*q]*c[(q-1, 2*q)] + normk2*c[(q-1,j-1)] 
        end
    end

    a2q = [c[(q-1, 2*q)] for q in 1:q_max]
    return a2q
end

@doc raw"""
    far_field_expansion_coeffs(stencil::Stencil{3, T}, m, PR::MPolyRing{T}) where T

Calculates an ``m``-term far-field expansion of the LGF using an algorithm from Martinnson and Rodin.
The result is an array `g` of multivariate polynomials in ``n`` so that
```math
G(n) = \frac{1}{4\pi |\bm{n}|} + \frac{1}{\pi}\sum_{k=1}^m \frac{g_k(\bm{n})}{|\bm{n}|^{6k+1}} + \mathcal{O}\left( |\bm{n}|^{-(2m + 3)} \right)
```

Note that once ``m`` exceeds 5 or 6, each additional term will roughly double the runtime of this function.

# Arguments
 - `stencil` : the 3D finite difference stencil
 - `m` : number of correction terms returned
 - `PR` : a multivariate `PolynomialRing` that defines the type of the output
"""
function far_field_expansion_coeffs(stencil::Stencil{3, T}, m, PR::MPolyRing = first(DefaultMultivariatePolyRing(stencil))) where {N, T}
    x = gens(PR)
    normx2 = sum(x.^2)
    a2q = martinsson_rodin_a2q(stencil, m, PR)

    # Now some combinatorics to evaluate g[q] = a[2q](∂) |x|^(2q-1) / fac
    # Mixed partials are not supported by AbstractAlgebra.jl, so we rely on the following.
    # For any multivariate polynomial p(x),
    # ∂ᵢ[p(x) |x|^n] = [∂ᵢp(x)|x|^2 + nxᵢp(x)] |x|^(n-2)
    g = zeros(PR, length(a2q))
    for (q, a) in enumerate(a2q)
        for (cff, ev) in zip(AbstractAlgebra.coefficients(a), exponent_vectors(a))
            p = PR(cff)
            n = 2*q-1
            for (i, xi) in enumerate(x)
                for deriv = 1:ev[i]
                    p = normx2*derivative(p, i) + n*xi*p
                    n -= 2
                end
            end
            g[q] += p
        end

        factor = (-1)^q * 2^(q+2) * fac(q) * semifactorial(2*q-1)
        g[q] = g[q] * (1 // factor)
    end

    return g
end

@doc raw"""
    far_field_expansion_coeffs(stencil::Stencil{2, T}, m, PR::MPolyRing{T}) where T

Calculates an ``m``-term far-field expansion of the 2D LGF using an algorithm from Martinnson and Rodin.
The result is an array `g` of multivariate polynomials in ``n`` so that
```math
G(\bm{n}) = -\frac{1}{2\pi}\left(\log |\bm{n}| + \gamma + \frac{\log 8}{2}\right) + \frac{1}{\pi}\sum_{k=1}^m \frac{g_k(\bm{n})}{|\bm{n}|^{6j}} + \mathcal{O}\left( |\bm{n}|^{-(m+1)} \right),
```
Note that once ``m`` exceeds 5 or 6, each additional term will roughly double the runtime of this function.

# Arguments
 - `stencil` : the 2D finite difference stencil
 - `m` : number of correction terms returned
 - `PR` : a multivariate `PolynomialRing` that defines the type of the output
"""
function far_field_expansion_coeffs(stencil::Stencil{2, T}, m, PR::MPolyRing = first(DefaultMultivariatePolyRing(stencil))) where T
    FF = FractionField(PR)
    x = gens(PR)
    normx2 = sum(x.^2)
    a2q = martinsson_rodin_a2q(stencil, m, PR)

    # Now some combinatorics to evaluate g[q] = a[2q](∂) (|x|^2q log |x|) / fac
    # Mixed partials are not supported by AbstractAlgebra.jl, so we rely on the following.
    # For any multivariate polynomial u(x), rational function v(x),
    # ∂ᵢ[u(x) log |x| + v(x)] = ∂ᵢu(x) log |x| + (xᵢ u(x)/|x|^2 + ∂ᵢv(x))
    g = zeros(PR, length(a2q))
    for (q, a) in enumerate(a2q)
        gq = FF(0)
        for (cff, ev) in zip(AbstractAlgebra.coefficients(a), exponent_vectors(a))

            u = cff * normx2^q
            v = FF(0)
            for (dim, xdim) in enumerate(x)
                for drvs = 1:ev[dim]
                    v = u * xdim // normx2 + derivative(v, dim)
                    u = derivative(u, dim)
                end
            end
            gq += v
        end

        factor = -(-1)^q * 2^(2*q+1) * fac(q)^2 # missing pi
        numer = (1 // factor) * gq * normx2^(3*q)
        @assert denominator(numer) == 1
        g[q] = numerator(numer)
    end

    return g
end
