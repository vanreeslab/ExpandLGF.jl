@doc raw"""
    tmax = inner_integral_errors(b, nmax; rtol=1e-16, atol=1e-16)

Calculate the region of validity for an inner integral expansion.

Specifically, when ``n \le n_{\mathrm{max}}`` and ``t \ge t_{\mathrm{max}}``, the last term in the 
expansion is smaller than both `atol` and `rtol * (first term)`.

# Arguments
 - `b` : vector of expansion coefficients (polynomials in n)
 - `nmax` : the maximal value of `n` that will be used in the expanions
 - `rtol` : relative tolerance on the last term in the series (pass 0 to ignore)
 - `atol` : absolute tolerance on the last term in the series (pass 0 to ignore)
"""
function inner_integral_errors(b, nmax; rtol=1e-16, atol=1e-16)
    @assert (rtol >= 0) && (atol >= 0) && (rtol + atol != 0) "Invalid rtol/atol combo."
    m = length(b)
    bm = abs(Float64(evaluate(b[m], nmax)))

    rtol_thresh = rtol == 0 ? -Inf : (bm / rtol)^(1/m)
    atol_thresh = atol == 0 ? -Inf : (bm / (atol * sqrt(4*pi)))^(2/(2*m+1))
    return max(rtol_thresh, atol_thresh)
end

@doc raw"""
    Tmax = outer_integral_errors(g, nmax; rtol=1e-16, atol=1e-16)

Calculate the region of validity for an outer integral expansion. The dimension (2D or 3D) is inferred from the type of `g`.

Specifically, when ``\|n\|_\infty ≤ n_{\mathrm{max}}`` and ``T \ge T_{\mathrm{max}}``, the last term in the 
expansion is smaller than both `atol` and `rtol * (first term)`.

# Arguments
 - `g` : vector of expansion coefficients (polynomials in n1, n2, n3)
 - `nmax` : the maximal value of the 3-vector `n` that will be used in the expanions
 - `rtol` : relative tolerance on the last term in the series (pass 0 to ignore)
 - `atol` : absolute tolerance on the last term in the series (pass 0 to ignore)
"""
function outer_integral_errors(g, nmax; rtol=1e-16, atol=1e-16)
    @assert (rtol >= 0) && (atol >= 0) && (rtol + atol != 0) "Invalid rtol/atol combo."
    m = length(g)
    dims = nvars(parent(g[1]))
    gm = abs(Float64(evaluate(g[m], nmax*ones(Int64, dims))))

    if (dims == 2)
        rtol_thresh = rtol == 0 ? -Inf : (gm / (rtol*g1))^(1/(m-1))
        atol_thresh = atol == 0 ? -Inf : (gm / (4*π*atol))^(1/m)
        return max(rtol_thresh, atol_thresh)
    else
        rtol_thresh = rtol == 0 ? -Inf : (gm / rtol)^(1/m)
        atol_thresh = atol == 0 ? -Inf : (gm / (atol * sqrt(16*π^3)))^(2/(2*m+1))
        return max(rtol_thresh, atol_thresh)
    end
end

@doc raw"""
    nmin = far_field_errors(g; rtol=1e-16, atol=1e-16, search=100)

Calculate the region of validity for a far-field expansion. The dimension (2D or 3D) is inferred from the type of `g`.

Specifically, when ``\|n\|_2 \ge n_{\mathrm{min}}``, the last term in the expansion is smaller than both 
`atol` and `rtol * (first term)`.

# Arguments
 - `g` : vector of far-field expansion coefficients (polynomials in n1, n2, n3)
 - `rtol` : relative tolerance on the last term in the series (pass 0 to ignore)
 - `atol` : absolute tolerance on the last term in the series (pass 0 to ignore)
 - `search` : finding `nmin` requires a bisection search over `n`. The `search` arg provides an upper limit.
"""
function far_field_errors(g; rtol=1e-16, atol=1e-16, search=100)
    dimension = nvars(parent(g[1]))
    @assert (rtol >= 0) && (atol >= 0) && (rtol + atol != 0) "Invalid rtol/atol combo."

    rtol = (rtol == 0) ? Inf : rtol
    atol = (atol == 0) ? Inf : atol
    m = length(g)

    # first_term(n) -> value of the first term in the expansions at x = (n, 0...)
    # last_term(n) -> value of the last term in the expansion at x = (n, 0...)
    first_term(n) = (dimension == 3) ? 
        1/(4*π*n) : 
        1/2π * (log(n) + MathConstants.γ + log(8)/2)

    last_term(n) = (dimension == 3) ? 
        abs(Float64(g[m](n, 0, 0))) / (π*Float64(n)^(6*m+1)) :
        abs(Float64(g[m](n, 0))) / (π*Float64(n)^(6*m)) 

    # bisection search
    nmin = 1
    nmax = search
    gmin = last_term(nmin)
    gmax = last_term(nmax)

    if (!(gmax < atol && gmax/first_term(nmax) < rtol))
        println("Search region error in far_field_errors()")
        return -1
    end

    while nmax - nmin > 1
        nnew = (nmin + nmax) ÷ 2
        gnew = last_term(nnew)

        if (gnew < atol && gnew/first_term(nnew) < rtol)
            nmax = nnew
            gmax = gnew
        else
            nmin = nnew
            gmin = gnew
        end
    end

    return nmax
end