@doc raw"""
    inner_integral_expansion(stencil::SplitStencil, terms)

Return a function which evaluates $I_{\sigma,n}(t)$ for large $t$ via a series expansion.

# Arguments
 - `stencil` : coefficients of the FD scheme
 - `terms` : the number of terms in the expansion

# Return
 - `out_fun` : a scalar function with signature `val = out_fun(n, t)`
"""
function inner_integral_expansion(stencil::SplitStencil{N, T}, terms) where {N, T}
    PR, _ = DefaultUnivariatePolyRing(stencil, "n")
    b = inner_integral_expansion_coeffs(stencil, terms, PR)
    return inner_integral_expansion(b)
end

@doc raw"""
    inner_integral_expansion(b)

Return a function which evaluates $I_{\sigma,n}(t)$ for large $t$ via a series expansion.

# Arguments
 - `b` : expansion coefficients, returned by `inner_integral_expansion_coeffs`

# Return
 - `out_fun` : a scalar function with signature `val = out_fun(n, t)`
"""
function inner_integral_expansion(b)
    @variables n, u
    expr = 1/sqrt(4*π) * (u + sum(poly_to_symbolic(bj, n)*u^(2*j+1) for (j, bj) in enumerate(b)))
    tmpfnc = build_function(expr, [n, u]; expression=Val{false})
    return (n, t) -> tmpfnc([n, 1/sqrt(t)])[1]
end

@doc raw"""
    outer_integral_expansion(stencil::SplitStencil{N, T}, terms) where {N, T}

Return a Julia function which evaluates the outer integral over $[t,\,\infty]$ for large $t$ via a series expansion.

# Arguments
 - `stencil` : coefficients of the FD scheme
 - `terms` : the number of terms in the expansion

# Return
 - `out_fun` : a scalar function with signature `val = out_fun(n, t)` for an N-Vector `n` and scalar `t`.
"""
function outer_integral_expansion(stencil::SplitStencil{N, T}, terms) where {N, T}
    MPR, _ = DefaultMultivariatePolyRing(stencil, "n")
    g = outer_integral_expansion_coeffs(stencil, terms, MPR)
    return (N == 2) ? outer_integral_2d_expansion(g) : outer_integral_3d_expansion(g)
end

@doc raw"""
    outer_integral_3d_expansion(g)

Return a Julia function which evaluates the 3D outer integral over $[T,\,\infty]$ for large $T$ via a series expansion.

# Arguments
 - `g` : expansion coefficients, the result of `outer_integral_3d_expansion_coeffs()`

# Return
 - `out_fun` : a scalar function with signature `val = out_fun([n1, n2, n3], T)`
"""
function outer_integral_3d_expansion(g)
    nsym = @variables n0, n1, n2
    @variables u
    expr = 1/sqrt(16*π^3) * (u + sum(poly_to_symbolic(gl, nsym)*u^(2*l+1) for (l, gl) in enumerate(g)))
    tmpfnc = build_function(expr, [n0, n1, n2, u]; expression=Val{false})
    return (nvec, T) -> tmpfnc([nvec[1], nvec[2], nvec[3], 1/sqrt(T)])[1]
end

@doc raw"""
    outer_integral_2d_expansion(g)

Return a Julia function which evaluates the 2D outer integral over $[T,\,\infty]$ for large $T$ via a series expansion.

# Arguments
 - `g` : expansion coefficients, the result of `outer_integral_2d_expansion_coeffs()`

# Return
 - `out_fun` : a scalar function with signature `val = out_fun([n1, n2], T)`
"""
function outer_integral_2d_expansion(g)
    nsym = @variables n1, n2
    @variables t
    expr = 1/(4*π) * sum(poly_to_symbolic(gl, nsym)*t^(-l) for (l, gl) in enumerate(g))
    tmpfnc = build_function(expr, [n1, n2, t]; expression=Val{false})
    return (nvec, T) -> tmpfnc([nvec[1], nvec[2], T])[1]
end 

"""
    far_field_expansion(stencil::Stencil{N, T}, nterms, target = Symbolics.JuliaTarget()) where {N, T}

Return a function which evaluates a far-field expansion of the LGF.

When `target = Symbolics.JuliaTarget()` (default), outputs a function `G = out_fun(n)` for scalar `G`, N-vector `n`.
When `target = Symbolics.CTarget()`, outputs a string containing a C function that evaluates the expansion.

# Arguments
 - `stencil` : coefficients of the finite difference stencil
 - `nterms` : number of terms in the expansion
 - `target` : choose either Symbolics.JuliaTarget()` or Symbolics.CTarget()
"""
 function far_field_expansion(stencil::Stencil{N, T}, nterms, target = Symbolics.JuliaTarget()) where {N, T}
    PR, _ = DefaultMultivariatePolyRing(stencil, "n")
    g = far_field_expansion_coeffs(stencil, nterms, PR)
    return (N == 2) ? far_field_2d_expansion(g, target) : far_field_3d_expansion(g, target)
end

"""
    far_field_3d_expression(g)

Return a Symbolics.jl expression for the expansion with coefficients in `g`.
"""
function far_field_3d_expression(g)
    nsym = @variables n1, n2, n3
    @variables normn
    expr = 1/(4*π*normn)
    for (q, gq) in enumerate(g)
        expr += poly_to_symbolic(gq, nsym) / (π * normn^(6*q + 1))
    end
    return n1, n2, n3, normn, expr
end

"""
    far_field_3d_expansion(g, target::JuliaTarget)

Return a Julia function which evaluates a far-field expansion of the 3D LGF.
The output has signature `G = out_fun(n)` for scalar `G`, 3-vector `n`
    
# Arguments
 - `g` : expansion coefficients, returned by `far_field_3d_expansion_coeffs`
"""
function far_field_3d_expansion(g, target::Symbolics.JuliaTarget)
    n1, n2, n3, normn, expr = far_field_3d_expression(g)
    tmpfnc, ~ = build_function([expr], [n1, n2, n3, normn]; expression=Val{false}, target = target)
    return (n) -> tmpfnc([n[1], n[2], n[3], sqrt(n[1]^2 + n[2]^2 + n[3]^2)])[1]
end

"""
    far_field_3d_expansion(g, target::CTarget)

Return a string containing C code that evaluates a far-field expansion of the 3D LGF.

# Arguments
 - `g` : far-field expansion coefficients, returned by `far_field_3d_expansion_coeffs`
"""
function far_field_3d_expansion(g, target::Symbolics.CTarget)
    n1, n2, n3, normn, expr = far_field_3d_expression(g)
    c_str = build_function([expr], n1, n2, n3, normn; expression=Val{true}, target = target,
        fname = :lgf, lhsname = :G, rhsnames = [:n1, :n2, :n3, :n])
    return replace(c_str, "π" => "M_PI");
end

"""
    far_field_2d_expression(g)

Return a Symbolics.jl expression for the expansion with coefficients in `g`.
"""
function far_field_2d_expression(g)
    nsym = @variables n1, n2
    @variables normn
    expr = -1/2π * (log(normn) + MathConstants.γ + log(8)/2)
    for (q, gq) in enumerate(g)
        expr += poly_to_symbolic(gq, nsym) / (π * normn^(6*q))
    end
    return n1, n2, normn, expr
end

"""
    far_field_2d_expansion(g, target::Symbolics.JuliaTarget)

Return a Julia function which evaluates a far-field expansion of the 2D LGF with coefficients in `g`.
The output has signature `G = out_fun(n)` for scalar `G`, 2-vector `n`

# Arguments
 - `g` : expansion coefficients, returned by `far_field_2d_expansion_coeffs`
"""
function far_field_2d_expansion(g, target::Symbolics.JuliaTarget)
    n1, n2, normn, expr = far_field_2d_expression(g)
    tmpfnc, ~ = build_function([expr], [n1, n2, normn]; expression=Val{false}, target = target)
    return (n) -> tmpfnc([n[1], n[2], sqrt(n[1]^2 + n[2]^2)])[1]
end

"""
    far_field_2d_expansion(g, target::Symbolics.CTarget)

Return a string containing C code that evaluates a far-field expansion of the 2D LGF.

# Arguments
 - `g` : expansion coefficients, returned by `far_field_2d_expansion_coeffs`
"""
function far_field_2d_expansion(g, target::Symbolics.CTarget)
    n1, n2, normn, expr = far_field_2d_expression(g)
    c_str = build_function([expr], n1, n2, normn; expression=Val{true}, target = target,
        fname = :lgf, lhsname = :G, rhsnames = [:n1, :n2, :n])
    return replace(c_str, "π" => "M_PI");
end


