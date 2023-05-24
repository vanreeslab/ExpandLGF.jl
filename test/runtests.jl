include("test_utils.jl")

# ---------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------
near_terms = 12
far_terms = 8
rtol = 0
atol = 1e-15
test_tol = 5e-14

for order in [2, 4, 6]
    println("--------------------------------------------------------")
    println("Testing 2D NearField for split stencil of order $(order)")
    println("--------------------------------------------------------")
    stencil = ExpandLGF.StandardDifference2D(order)
    nmax = 16
    NF = ExpandLGF.NearField(stencil, nmax, near_terms, rtol, atol)

    # Integrate to get indices i1 >= i2
    @printf("Testing 2D NearField for indices [0, %d]^2 (may take up to a minute)\n", nmax)
    N = nmax
    lgf = OffsetArray(zeros(Float64, (N+1, N+1)), 0:N, 0:N)
    @time for i1 in 0:N, i2 in 0:i1
        lgf[i1, i2] = NF([i1, i2])
    end

    # Test the max residual within the block
    test = generate_full_test_block(stencil, lgf)
    test_lgf_block(stencil, test, test_tol)
end

for order in [2, 4, 6]
    println("--------------------------------------------------------")
    println("Testing 2D FarField for split stencil of order $(order)")
    println("--------------------------------------------------------")
    stencil = ExpandLGF.StandardDifference2D(order)
    FF, nmin = ExpandLGF.FarField(stencil, far_terms, rtol, atol)

    # Integrate to get indices i1 >= i2
    org = [nmin, nmin]
    N = 12
    rngs = (org[1] .+ (0:N), org[2] .+ (0:N))
    lgf = OffsetArray(zeros(Float64, (N+1, N+1)), rngs...)
    @printf("Testing 2D FarField for a block of size %d^3 with origin [%d, %d]\n", N+1, org[1], org[2])
    @time for i1 in rngs[1], i2 in rngs[2]
        lgf[i1, i2] = FF([i1, i2])
    end

    test_lgf_block(stencil, lgf, test_tol)
end

for order in [2, 4, 6]
    println("--------------------------------------------------------")
    println("Testing 2D ExpandLGF for split stencil of order $(order)")
    println("--------------------------------------------------------")
    stencil = ExpandLGF.StandardDifference2D(order)
    ELGF = ExpandLGF.EvaluateLGF(stencil, near_terms, far_terms, rtol, atol)

    # Integrate to get indices i1 >= i2
    N = 32
    lgf = OffsetArray(zeros(Float64, (N+1, N+1)), 0:N, 0:N)
    @printf("Testing 2D EvaluateLGF for indices [0, %d]^2 (may take up to a minute)\n", N)
    @time for i1 in 0:N, i2 in 0:i1
        lgf[i1, i2] = ELGF([i1, i2])
    end

    # Test the max residual within the block
    test = generate_full_test_block(stencil, lgf)
    test_lgf_block(stencil, test, test_tol)
end

for order in [2, 4, 6]
    println("--------------------------------------------------------")
    println("Testing 3D NearField for split stencil of order $(order)")
    println("--------------------------------------------------------")
    stencil = ExpandLGF.StandardDifference3D(order)
    nmax = 16
    NF = ExpandLGF.NearField(stencil, nmax, near_terms, rtol, atol)

    # Integrate to get indices i1 >= i2 >= i3
    @printf("Testing 3D NearField for indices [0, %d]^3 (may take up to a minute)\n", nmax)
    N = nmax
    lgf = OffsetArray(zeros((N+1, N+1, N+1)), 0:N, 0:N, 0:N)
    @time for i1 in 0:N, i2 in 0:i1, i3 in 0:i2 
        lgf[i1, i2, i3] = NF([i1, i2, i3])
    end

    # Test the max residual within the block
    test = generate_full_test_block(stencil, lgf)
    test_lgf_block(stencil, test, test_tol)
end

for (stencil, order) in [(ExpandLGF.LeftMehrstellen4(), 4), (ExpandLGF.LeftMehrstellen6(), 6)]
    println("--------------------------------------------------------")
    println("Testing 3D NearField for left Mehrstellen stencil of order $(order)")
    println("--------------------------------------------------------")
    mehrstellen_atol = 1e-8
    mehrstellen_evals = 1e8
    nmax = 4
    NF = ExpandLGF.NearField(stencil, 0, mehrstellen_atol, mehrstellen_evals)

    @printf("Testing 3D NearField for indices [0, %d]^3 (may take up to a minute)\n", nmax)
    N = nmax
    lgf = OffsetArray(zeros((N+1, N+1, N+1)), 0:N, 0:N, 0:N)
    @time for i1 in 0:N, i2 in 0:i1, i3 in 0:i2 
        lgf[i1, i2, i3] = NF([i1, i2, i3])
    end

    test = generate_full_test_block(stencil, lgf)
    test_lgf_block(stencil, test, mehrstellen_atol)
end

for (stencil, order) in [(ExpandLGF.Mehrstellen4(), 4), (ExpandLGF.Mehrstellen6(), 6)]
    println("--------------------------------------------------------")
    println("Testing 3D NearField for Mehrstellen stencil of order $(order)")
    println("--------------------------------------------------------")
    mehrstellen_atol = 1e-8
    mehrstellen_rtol = 0
    mehrstellen_evals = 1e8
    nmax = 4
    NF = ExpandLGF.NearField(stencil, mehrstellen_rtol, mehrstellen_atol, mehrstellen_evals)

    @printf("Testing 3D NearField for indices [0, %d]^3 (may take up to a minute)\n", nmax)
    N = nmax
    lgf = OffsetArray(zeros((N+1, N+1, N+1)), 0:N, 0:N, 0:N)
    @time for i1 in 0:N, i2 in 0:i1, i3 in 0:i2 
        lgf[i1, i2, i3] = NF([i1, i2, i3])
    end

    test = generate_full_test_block(ExpandLGF.left_stencil(stencil), lgf)
    test_lgf_block(stencil, test, mehrstellen_atol)
end

for order in [2, 4, 6]
    println("--------------------------------------------------------")
    println("Testing 3D FarField for split stencil of order $(order)")
    println("--------------------------------------------------------")
    stencil = ExpandLGF.StandardDifference3D(order)
    FF, nmin = ExpandLGF.FarField(stencil, far_terms, rtol, atol)

    org = [nmin, nmin, nmin]
    N = 12
    rngs = (org[1] .+ (0:N), org[2] .+ (0:N), org[3] .+ (0:N))
    lgf = OffsetArray(zeros(Float64, (N+1, N+1, N+1)), rngs...)
    @printf("Testing 3D FarField for a block of size %d^3 with origin [%d, %d, %d]\n", N+1, org[1], org[2], org[3])
    @time for i1 in rngs[1], i2 in rngs[2], i3 in rngs[3] 
        lgf[i1, i2, i3] = FF([i1, i2, i3])
    end

    test_lgf_block(stencil, lgf, test_tol)
end

for (stencil, order) in [(ExpandLGF.LeftMehrstellen4(), 4), (ExpandLGF.LeftMehrstellen6(), 6)]
    println("--------------------------------------------------------")
    println("Testing 3D FarField for left Mehrstellen stencil of order $(order)")
    println("--------------------------------------------------------")
    FF, nmin = ExpandLGF.FarField(stencil, far_terms, rtol, atol)

    org = [nmin, nmin, nmin]
    N = 12
    rngs = (org[1] .+ (0:N), org[2] .+ (0:N), org[3] .+ (0:N))
    lgf = OffsetArray(zeros(Float64, (N+1, N+1, N+1)), rngs...)
    @printf("Testing 3D FarField for a block of size %d^3 with origin [%d, %d, %d]\n", N+1, org[1], org[2], org[3])
    @time for i1 in rngs[1], i2 in rngs[2], i3 in rngs[3] 
        lgf[i1, i2, i3] = FF([i1, i2, i3])
    end

    test_lgf_block(stencil, lgf, test_tol)
end

for (stencil, order) in [(ExpandLGF.Mehrstellen4(), 4), (ExpandLGF.Mehrstellen6(), 6)]
    println("--------------------------------------------------------")
    println("Testing 3D FarField for combined Mehrstellen stencil of order $(order)")
    println("--------------------------------------------------------")
    FF, nmin = ExpandLGF.FarField(stencil, far_terms, rtol, atol)

    org = [nmin, nmin, nmin]
    N = 12
    rngs = (org[1] .+ (0:N), org[2] .+ (0:N), org[3] .+ (0:N))
    lgf = OffsetArray(zeros(Float64, (N+1, N+1, N+1)), rngs...)
    @printf("Testing 3D FarField for a block of size %d^3 with origin [%d, %d, %d]\n", N+1, org[1], org[2], org[3])
    @time for i1 in rngs[1], i2 in rngs[2], i3 in rngs[3] 
        lgf[i1, i2, i3] = FF([i1, i2, i3])
    end

    test_lgf_block(stencil, lgf, test_tol)
end

for order in [2, 4, 6]
    println("--------------------------------------------------------")
    println("Testing 3D ExpandLGF for split stencil of order $(order)")
    println("--------------------------------------------------------")
    stencil = ExpandLGF.StandardDifference3D(order)
    ELGF = ExpandLGF.EvaluateLGF(stencil, near_terms, far_terms, rtol, atol)

    # Integrate to get indices i1 >= i2 >= i3
    N = 32
    lgf = OffsetArray(zeros(Float64, (N+1, N+1, N+1)), 0:N, 0:N, 0:N)
    @printf("Testing 3D EvaluateLGF for indices [0, %d]^3 (may take up to a minute)\n", N)
    @time for i1 in 0:N, i2 in 0:i1, i3 in 0:i2 
        lgf[i1, i2, i3] = ELGF([i1, i2, i3])
    end

    # Test the max residual within the block
    test = generate_full_test_block(stencil, lgf)
    test_lgf_block(stencil, test, test_tol)
end