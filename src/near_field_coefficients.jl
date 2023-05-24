@doc raw"""
    inner_integral_expansion_coeffs(stencil::SplitStencil, n, m)

Return an ``m``-term expansion of the inner integral for fixed ``n`` and large ``t``.
The results is an array `b` of univariate polynomials in `n` so that 

```math
I_{\sigma, n}(t) = \frac{1}{\sqrt{4\pi t}} \sum_{k=0}^m b_k(n) t^{-k} + \mathcal{O}\left( t^{-\frac{2m+3}{2}} \right)
```

# Arguments
 - `stencil` : the finite difference stencil
 - `m` : the number of terms in the expansion
 - `PR` : a univariate `PolynomialRing`
"""    
function inner_integral_expansion_coeffs(stencil::SplitStencil, m, PR::PolyRing = first(DefaultUnivariatePolyRing(stencil)))

    # power series ring in k, with coeffs in QQ[n, t]
    MPR, (nm, t) = PolynomialRing(QQ, ["n", "t"])
    UPSR, k = PowerSeriesRing(MPR, 4*m+1, "k"; model=:capped_absolute)

    # get expansion of exp(-t(σ - k²))*cos(nk)
    σ = sum((∂σ(coefficients(stencil), i) // fac(i) * k^i for i = 4:2:4*m)) 
    cosnk = sum(((-1)^i // fac(2i) * (nm*k)^(2i)  for i = 0:2*m))
    int = exp(-t*σ) * cosnk # int = sum_j a[2j]/(2j)! k^(2j)
    
    # integrate each term wrt k and recombine
    nu = gen(PR)
    b = zeros(PR, 2*m)
    for j in 1:2*m
        poly = coeff(int, 2j) * semifactorial(2j - 1) // 2^j # = a2j(n,t)/((2j)!! 2^j) 
        for it in 1:length(poly)
            en, et = exponent_vector(poly, it)
            b[j - et] += coeff(poly, it)*nu^en
        end
    end

    # return the fully resolved b[2j](n), half of those present 
    return b[1:m]
end

@doc raw"""
    outer_integral_expansion_coeffs(stencil::SplitStencil{3, T}, n, m) where T

Return an ``m``-term expansion of the 3D outer integral for large ``T``.
The results is an array ``g`` of mulitvariate polynomials in `n` so that 

```math
\int_T^\infty I_{\sigma, n_1}(t) I_{\sigma, n_2}(t) I_{\sigma, n_3}(t) \,\mathrm{d}t \approx \frac{1}{\sqrt{16\pi^3T}} \sum_{k=0}^m g_k(\bm{n}) T^{-k} + \mathcal{O}\left( T^{-\frac{2m+3}{2}} \right)
```

# Arguments
 - `stencil` : the finite difference stencil
 - `n` : a vector of 3 generators from a multivariate `PolynomialRing`
 - `m` : the number of terms in the expansion
"""
function outer_integral_expansion_coeffs(stencil::SplitStencil{3, T}, nterms::Int, MPR::MPolyRing) where T
    UPR, _ = DefaultUnivariatePolyRing(stencil, "n")
    b = inner_integral_expansion_coeffs(stencil, nterms, UPR)
    return outer_integral_3d_expansion_coeffs(b, MPR)
end

"""
    outer_integral_3d_expansion_coeffs(b, n)

Return a power series expansion of the 3D outer integral based on an expansion of the inner integral.

# Arguments
 - `b` : inner expansion, the result of `inner_integral_expansion_coeffs()`
 - `PR` : a multivariate `PolynomialRing` with 3 generators
"""
function outer_integral_3d_expansion_coeffs(b, PR::MPolyRing)
    
    n1, n2, n3 = gens(PR)
    m = length(b)

    b1 = map((bj) -> AbstractAlgebra.subst(bj, n1), b)
    b2 = map((bj) -> AbstractAlgebra.subst(bj, n2), b)
    b3 = map((bj) -> AbstractAlgebra.subst(bj, n3), b)
    g = zeros(PR, m)
    for j1 = 0:m, j2 = 0:m, j3 = 0:m
        l = j1 + j2 + j3
        if l > m || l == 0
            continue
        end
        g[l] += (j1 == 0 ? 1 : b1[j1]) *
                (j2 == 0 ? 1 : b2[j2]) *  
                (j3 == 0 ? 1 : b3[j3]) * 
                (1 // (2l + 1))
    end
    return g
end

@doc raw"""
    outer_integral_expansion_coeffs(stencil::SplitStencil{2, T}, n, m) where T

Return an ``m``-term expansion of the 2D outer integral for large ``T``.
The results is an array ``g`` of mulitvariate polynomials in `n` so that 

```math
\int_T^\infty [I_{\sigma, n_1}(t) I_{\sigma, n_2}(t) - I_{\sigma, 0}(t)^2] \,\mathrm{d}t \approx \frac{1}{4\pi} \sum_{k=1}^m g_k(\bm{n}) T^{-k} + \mathcal{O}\left( T^{-(m+1)} \right)
```

# Arguments
 - `stencil` : the finite difference stencil
 - `n` : a vector of 2 generators from a multivariate `PolynomialRing`
 - `m` : the number of terms in the expansion
"""
function outer_integral_expansion_coeffs(stencil::SplitStencil{2, T}, nterms::Int, MPR::MPolyRing) where T
    UPR, _ = DefaultUnivariatePolyRing(stencil, "n")
    b = inner_integral_expansion_coeffs(stencil, nterms, UPR)
    return outer_integral_2d_expansion_coeffs(b, MPR)
end

"""
    outer_integral_2d_expansion_coeffs(b, PR::MPolyRing)

Return a power series expansion of the 2D outer integral based on an expansion of the inner integral.

# Arguments
 - `b` : inner expansion, the result of `inner_integral_expansion_coeffs()`
 - `PR` : a multivariate `PolynomialRing` in two variables
"""
function outer_integral_2d_expansion_coeffs(b, PR::MPolyRing)
    
    n1, n2 = gens(PR)
    m = length(b)

    b0 = map((bj) -> AbstractAlgebra.subst(bj, 0), b)
    b1 = map((bj) -> AbstractAlgebra.subst(bj, n1), b)
    b2 = map((bj) -> AbstractAlgebra.subst(bj, n2), b)

    g = zeros(PR, m)
    for j1 = 0:m, j2 = 0:m
        l = j1 + j2
        if l > m || l == 0
            continue
        end
        g[l] += (1 // l) * 
                ((j1 == 0 ? 1 : b1[j1]) * (j2 == 0 ? 1 : b2[j2]) -
                 (j1 == 0 ? 1 : b0[j1]) * (j2 == 0 ? 1 : b0[j2]))
    end
    return g
end