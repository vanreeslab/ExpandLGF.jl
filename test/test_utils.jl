using ExpandLGF
using ExpandLGF: sorted_indices_2D, sorted_indices_3D
using OffsetArrays
using Test
using Printf

"Return a zero OffsetArray obtained by adding `stencil_width` layers of points to a `block`"
function allocate_resized_array(block::OffsetArray, expand_low::Int, expand_high::Int = expand_low)
    block_is_resizable = all(block_width + expand_low + expand_high > 0 for block_width in size(block))
    @assert(block_is_resizable, "cannot resize a block of size $(size(block)) by ($(expand_low), $(expand_high)) points")
    
    expand_range = (range) -> (first(range) - expand_low):(last(range) + expand_high)
    output_ranges = map(expand_range, axes(block))
    output_size = (block_width + expand_low + expand_high for block_width in size(block)) 

    return OffsetArray(zeros(output_size...), output_ranges...)
end

"If `out` contains the origin, add 1.0 to that entry"
function handle_origin!(out)
    output_contains_origin = all(0 in axis for axis in axes(out))
    if output_contains_origin
        origin = repeat([0], ndims(out))
        out[origin...] += 1.0
    end
    return out
end

"Subtract the RHS stencil of a MehrstellenStencil from the appropraite indices of the array `out`"
function handle_origin!(out, stencil::ExpandLGF.MehrstellenStencil{N, T}) where {N, T}
    rs = ExpandLGF.right_stencil(stencil)
    for (index, value) in zip(ExpandLGF.indices(rs), ExpandLGF.values(rs))
        output_contains_index = all(index[i] in axes(out, i) for i in 1:N)
        if output_contains_index
            out[index...] += value
        end
    end
    return out
end

"""
    out = apply_stencil(stencil::ExpandLGF.SplitStencil, block)

Return the result of applying `stencil` to the entries of `block`.

# Arguments
 - `stencil` : any stencil type from the ExpandLGF package
 - `block` : an OffsetArray to apply the stencil too

 # Returns
  - `out` : an OffsetArray containing entries for points where `stencil` can be evaluated
"""
function apply_stencil(stencil::ExpandLGF.SplitStencil, block)

    # Decode stencil
    w = ExpandLGF.width(stencil)
    stencil_coefficients = OffsetArray(zeros(2*w + 1), -w:w)
    for (i, ci) in enumerate(ExpandLGF.coefficients(stencil))
        stencil_coefficients[-i] = ci
        stencil_coefficients[0] += -2*ci
        stencil_coefficients[+i] = ci
    end

    # Apply stencil
    out = allocate_resized_array(block, -w)
    for I in CartesianIndices(out)
        for dim = 1:ndims(block)
            for l = -w:w
                J = [Tuple(I)...]
                J[dim] += l
                out[I] += stencil_coefficients[l] * block[J...]
            end
        end
    end
    return out
end

function apply_stencil(stencil::ExpandLGF.FullStencil, block)
    
    # Decode stencil
    width = ExpandLGF.width(stencil)
    indices = ExpandLGF.indices(stencil)
    values = ExpandLGF.values(stencil)

    out = allocate_resized_array(block, -width)
    for I in CartesianIndices(out)
        for (J, s) in zip(indices, values)
            point = [Tuple(I)...] + [J...]
            out[I] += s * block[point...]
        end
    end
    return out
end

"Given LGF entries for 0 <= i1 <= ... <= iN, generate a block with -w <= i <= N in each dimension"
function generate_full_test_block(stencil, lgf)
    w = ExpandLGF.width(stencil)
    test = allocate_resized_array(lgf, w, 0)
    for I in CartesianIndices(test)
        J = sort(abs.([Tuple(I)...]), rev=true)
        test[I] = lgf[J...]
    end
    return test
end

"Assert that the max residual after applying `stencil` to `test` is below `tol`"
function test_lgf_block(stencil::ExpandLGF.Stencil, test, tol)
    println("Doing the standard version")
    residual = apply_stencil(stencil, test)
    handle_origin!(residual)
    max_residual, index = findmax(abs, residual)
    @printf("Max residual: %1.3e", max_residual)
    print(" at index $(Tuple(index))\n\n")
    @test max_residual < tol
end

"Assert that the max residual after applying `stencil` to `test` is below `tol`"
function test_lgf_block(stencil::ExpandLGF.MehrstellenStencil, test, tol)
    println("Doing the mehrstellen version")
    residual = apply_stencil(ExpandLGF.left_stencil(stencil), test)
    handle_origin!(residual, stencil)
    max_residual, index = findmax(abs, residual)
    @printf("Max residual: %1.3e", max_residual)
    print(" at index $(Tuple(index))\n\n")
    @test max_residual < tol
end

