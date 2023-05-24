#############
# NearField #
#############

"""
A callable structure that evaluates near-field values through quadrature.
Call with `(::NearFieldEvaluator)(n)` for an integer vector `n`. 
"""
abstract type NearFieldEvaluator end

"A callable structure that evaluates near-field values for a dimension-split stencil"
struct NearFieldSplit <: NearFieldEvaluator
    dimension::Int64
    stencil::SplitStencil
    rtol::Float64
    atol::Float64
    inner_expansion::Function
    outer_expansion::Function
    inner_thresh::Float64
    outer_thresh::Float64
end

Base.show(io::IO, NF::NearFieldSplit) = print(io, 
"""ExpandLGF.NearFieldSplit with
    dimension = $(NF.dimension)
    stencil = $(NF.stencil)
    rtol = $(NF.rtol)
    atol = $(NF.atol)
    inner_thresh = $(NF.inner_thresh)
    outer_thresh = $(NF.outer_thresh)
""")

@doc raw"""
    NearField(stencil::SplitStencil, nmax, terms, rtol, atol)

Return a callable struct which evaluates the LGF for a `SplitStencil` for ``\|n\|_\infty \le {n_\mathrm{max}}`` with prescribed error tolerance.
Uses a combination of direct quadrature and series expansions.

# Arguments
 - `stencil` : the finite difference stencil
 - `nmax` : the largest component of `n` that will be evaluated
 - `terms` : number of terms used in integral expansions. Affects performance.
 - `rtol`, `atol` : relative and absolute error tolerance for LGF values
"""
function NearField(stencil::SplitStencil{N, T}, nmax, terms, rtol, atol) where {N, T}

    # get expansion coeffs (1 extra term for error estimation)
    UPR, _ = DefaultUnivariatePolyRing(stencil)
    inner_coeffs = inner_integral_expansion_coeffs(stencil, terms+1, UPR)

    MPR, _ = DefaultMultivariatePolyRing(stencil)
    outer_coeffs = (N == 2) ? 
        outer_integral_2d_expansion_coeffs(inner_coeffs, MPR) :
        outer_integral_3d_expansion_coeffs(inner_coeffs, MPR)

    # based on the tolerances, get thresholds
    inner_thresh = inner_integral_errors(inner_coeffs, nmax; rtol = rtol, atol = atol)
    outer_thresh = outer_integral_errors(outer_coeffs, nmax; rtol = rtol, atol = atol)
    
    # build callable expansions
    inner_expansion = inner_integral_expansion(inner_coeffs[1:terms])
    outer_expansion = (N == 2) ? 
        outer_integral_2d_expansion(outer_coeffs[1:terms]) :
        outer_integral_3d_expansion(outer_coeffs[1:terms])

    # return callable struct
    return NearFieldSplit(N, stencil, rtol, atol, inner_expansion, outer_expansion, inner_thresh, outer_thresh)
end

function (NF::NearFieldSplit)(n)
    @assert length(n) == NF.dimension

    function inner_integral_(n, t)
        if t > NF.inner_thresh
            return NF.inner_expansion(n, t)
        else
            return inner_integral_quad(NF.stencil, n, t)
        end
    end

    outer_integrand_ = (NF.dimension == 2) ? 
        (t) -> (inner_integral_(n[1], t)*inner_integral_(n[2], t) - inner_integral_(0, t)^2) :
        (t) -> (prod(inner_integral_(ni, t) for ni in n))

    num, ~ = quadgk(outer_integrand_, 0, NF.outer_thresh; rtol = NF.rtol, atol = NF.atol, order=21)
    anl = NF.outer_expansion(n, NF.outer_thresh)
    
    return num + anl
end

"A callable structure that evaluates near-field values via quadrature for a full or Mehrstellen stencil"
struct NearFieldFull <: NearFieldEvaluator
    dimension::Int
    symbol::Function
    rtol::Float64
    atol::Float64
    maxevals::Int
end

Base.show(io::IO, NF::NearFieldFull) = print(io, 
"""ExpandLGF.NearFieldFull with
    dimension = $(NF.dimension)
    rtol = $(NF.rtol)
    atol = $(NF.atol)
    maxevals = $(NF.maxevals)
""")

@doc raw"""
    NearField(stencil::Stencil, rtol, atol, maxevals)

Return a callable struct which evaluates the LGF for any `Stencil` with prescribed error tolerance.
Uses direct quadrature.

# Arguments
 - `stencil` : the finite difference stencil
 - `nmax` : the largest component of `n` that will be evaluated
 - `rtol`, `atol` : relative and absolute error tolerance for LGF values
 - `maxevals` : the max number of symbol evaluations used to compute each LGF value
"""
function NearField(stencil::Stencil{N, T}, rtol=0, atol=1e-9, maxevals=10^8) where {N, T}
    return NearFieldFull(N, ExpandLGF.compiled_symbol(stencil), rtol, atol, maxevals)
end

function (NF::NearFieldFull)(n)
    if (NF.dimension == 2)
        f = (k) -> -1/π^2 * (cos(n[1]*k[1])*cos(n[2]*k[2]) - 1.0) / NF.symbol(k)
        integral, error = hcubature(f, [0.,0.], [π, π]; rtol=NF.rtol, atol=NF.atol, maxevals=NF.maxevals)
        return integral
    else
        f = (k) -> -1/π^3 * cos(n[1]*k[1])*cos(n[2]*k[2])*cos(n[3]*k[3]) / NF.symbol(k)
        integral, error = hcubature(f, [0.,0.,0.], [π, π, π]; rtol=NF.rtol, atol=NF.atol, maxevals=NF.maxevals)
        return integral
    end
end

############
# FarField #
############

"""
A callable struct for far field evaluation. Call via `(::FarFieldExpansion)(n)` for a 2D or 3D integer vector `n`. 
"""
struct FarFieldExpansion
    dimension::Int64
    expansion_::Function
end

Base.show(io::IO, FF::FarFieldExpansion) = print(io, "ExpandLGF.FarFieldExpansion for a $(FF.dimension)D stencil")

@doc raw"""
    far_field FarField(stencil::Stencil, nterms)

Evaluate the LGF for `stencil` via a far-field series expansion with `ntemrs` terms.
"""
function FarField(stencil::Stencil{N, T}, nterms) where {N, T}
    return FarFieldExpansion(N, far_field_expansion(stencil, nterms))
end

@doc raw"""
    far_field, nmin = FarField(stencil::Stencil, nterms, rtol, atol)

Evaluate the LGF for ``\|n\|_2 \ge n_{\mathrm{min}}`` via series expansion with prescribed error tolerance.
Returns a callable struct and the threshold `nmin` beyond which the expansion achieves the given tolerance.

# Arguments
 - `stencil` - the finite difference stencil
 - `dimensions` : either 2 or 3 to indicate a 2D or 3D LGF
 - `nterms` - number of terms in the expansion
 - `rtol`, `atol` - relative and absolute error tolerance on the evaluation (pass 0 to ignore either).
"""
function FarField(stencil::Stencil{N, T}, nterms, rtol, atol) where {N, T}
    MPR, _ = DefaultMultivariatePolyRing(stencil, "n")
    g = far_field_expansion_coeffs(stencil, nterms + 1, MPR)
    far_cutoff = far_field_errors(g; rtol = rtol, atol = atol)
    far_expansion = (N == 2) ?
        far_field_2d_expansion(g[1:nterms], Symbolics.JuliaTarget()) :
        far_field_3d_expansion(g[1:nterms], Symbolics.JuliaTarget())
    return FarFieldExpansion(N, far_expansion), far_cutoff
end

"""
    (FF::FarFieldExpansion)(n)

Evaluate an LGF via `FarField` expansion at a 2D or 3D integer vector `n`. 
"""
function(FF::FarFieldExpansion)(n)
    @assert length(n) == FF.dimension
    return FF.expansion_(n)
end

###############
# EvaluateLGF #
###############

"""
A callable struct wrapping near-field evaluation, far-field evaluation, and a cutoff between the two.
Call via `(::EvaluateLGF)(n)` for a 2D or 3D integer vector `n`. 
"""
struct EvaluateLGF
    dimension::Int64
    nf_::NearFieldEvaluator
    ff_::FarFieldExpansion
    rad_::Int64
    offset_::Float64
end

near_field(ELGF::EvaluateLGF) = ELGF.nf_
far_field(ELGF::EvaluateLGF) = ELGF.ff_
near_cutoff(ELGF::EvaluateLGF) = ELGF.rad_

function (ELGF::EvaluateLGF)(n)
    @assert length(n) == ELGF.dimension
    if sum(n.^2) >= ELGF.rad_^2
        return ELGF.ff_(n) + ELGF.offset_
    else
        return ELGF.nf_(n)
    end
end

"""
    EvaluateLGF(stencil::SplitStencil, near_terms, far_terms, rtol, atol)

Return a callable struct that evaluates the LGF by switching between near-field and far-field evaluation strategies.

# Arguments
 - `stencil` : the split finite difference stencil
 - `near_terms` : number of inner and outer integral expansion terms.
 - `far_terms` : number of far-field expansion terms.
 - `rtol`, `atol` : relative and absolute error tolerance for the evaluation.
"""
function EvaluateLGF(stencil::SplitStencil{N, T}, near_terms, far_terms, rtol, atol) where {N, T}

    far_field, far_cutoff = FarField(stencil, far_terms, rtol, atol)
    near_field = NearField(stencil, far_cutoff, near_terms, rtol, atol)
    offset = (N == 2) ? near_field([far_cutoff, 0]) - far_field([far_cutoff, 0]) : 0.0

    return EvaluateLGF(N, near_field, far_field, far_cutoff, offset)
end

"""
    EvaluateLGF(stencil::Union{FullStencil, MehrstellenStencil}, far_terms, rtol, atol, maxevals)

Return a callable struct that evaluates the LGF for a `FullStencil` or `MehrstellenStencil`, switching between near-field or far-field evaluation strategies.

# Arguments
 - `stencil` : the finite difference stencil
 - `far_terms` : number of far-field expansion terms.
 - `rtol`, `atol` : relative and absolute error tolerance for the evaluation.
 - `maxevals` : the max number of symbol evaluations used to compute near-field LGF values
"""
function EvaluateLGF(stencil::Union{FullStencil{N, T}, MehrstellenStencil{N, T}}, far_terms, rtol, atol, maxevals) where {N, T}

    far_field, far_cutoff = FarField(stencil, far_terms, rtol, atol)
    near_field = NearField(stencil, rtol, atol, maxevals)
    offset = (N == 2) ? near_field([far_cutoff, 0]) - far_field([far_cutoff, 0]) : 0.0

    return EvaluateLGF(N, near_field, far_field, far_cutoff, offset)
end

function Base.show(io::IO, ELGF::EvaluateLGF)
    if (ELGF.dimension == 3)
        print(io, "ExpandLGF.EvaluateLGF for a 3D stencil, cutoff = $(ELGF.rad_)")
    else
        print(io, "ExpandLGF.EvaluateLGF for a 2D stencil, cutoff = $(ELGF.rad_) and offset = $(ELGF.offset_)")
    end
end