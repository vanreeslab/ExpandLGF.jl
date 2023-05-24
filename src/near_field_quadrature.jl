"""
    inner_integrand(stencil::SplitStencil, n, t, k)

Return `exp(-tσ)*cos(nk)`, where σ is the symbol of the 1D finte difference
"""
function inner_integrand(stencil::SplitStencil, n, t, k)
    σ = 4*sum(cj*sin(j*k/2)^2 for (j, cj) in enumerate(coefficients(stencil)))
    return exp(-t*σ)*cos(n*k)
end

"""
    inner_integral_quad(cffs, n, t; rtol = 1e-15, atol = 1e-15, order = 21)

Evaluate `I[σ](n,t)` through numerical quadrature with the QuadGK package.

# Arguments
 - rtol : relative tolerance 
 - atol : absolute tolerance
 - order : quadrature order passed to QGK
"""
function inner_integral_quad(stencil::SplitStencil, n, t; rtol = 1e-15, atol = 1e-15, order=21)
    integrand_(k::Float64) = inner_integrand(stencil, n, t, k)
    int = quadgk(integrand_, 0, π; rtol = rtol, atol = π*atol, order=order)
    return int[1] / π
end

"""
    evaluate_symbol(stencil::FullStencil, k)

Evaluate the symbol `σ(k)` for a `FullStencil`.
"""
function evaluate_symbol(stencil::FullStencil, k)
    out = 0
    for (α, s) in zip(indices(stencil), values(stencil))
        out += s * prod(cos(abs(αi)*ki) for (αi, ki) in zip(α, k))
    end
    return out
end

"""
    polynomial_symbol(stencil::FullStencil, MPR::MPolyRing)

Return a polynomial `p` so that `σ(k) = p(sin^2(k1/2), sin^2(k2/2), ...)`,
where `σ(k)` is the symbol for a `FullStencil`.
"""
function polynomial_symbol(stencil::FullStencil, MPR::MPolyRing)
    # cache polynomials s.t. cos(n*k) = P[n+1](sin(k/2)^2) 
    UPR, x = DefaultUnivariatePolyRing(stencil)
    w = width(stencil)
    P = OffsetArray([AbstractAlgebra.chebyshev_t(n, 1 - 2*x) for n in 0:w], 0:w)
    # compute output
    k = gens(MPR)
    out = zero(MPR)
    for (α, s) in zip(indices(stencil), values(stencil))
        out += s * prod(P[abs(αi)](ki) for (αi, ki) in zip(α, k))
    end
    return out
end

"""
    compiled_symbol(stencil)

Return a compiled Julia function that evaluates the symbol for a `FullStencil` or `MehrstellenStencil`.
"""
function compiled_symbol(stencil::FullStencil)
    MPR, _ = DefaultMultivariatePolyRing(stencil)
    polyint = polynomial_symbol(stencil, MPR)

    @syms k1 k2 k3
    expr = poly_to_symbolic(polyint, [k1, k2, k3])
    compiled_function = build_function(expr, k1, k2, k3; expression=Val{false})
    return (k) -> compiled_function((sin.(k/2).^2)...)
end


function compiled_symbol(stencil::MehrstellenStencil)
    MPR, _ = DefaultMultivariatePolyRing(stencil)
    left_poly = polynomial_symbol(left_stencil(stencil), MPR)
    right_poly = polynomial_symbol(right_stencil(stencil), MPR)

    @syms k1 k2 k3
    expr = poly_to_symbolic(left_poly, [k1, k2, k3]) / poly_to_symbolic(right_poly, [k1, k2, k3])
    compiled_function = build_function(expr, k1, k2, k3; expression=Val{false})
    return (k) -> compiled_function((sin.(k/2).^2)...)
end

