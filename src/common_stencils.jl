"Return the stencil for the 4th order left-hand side Mehrstellen operator"
function LeftMehrstellen4()
    return FullStencil(Dict((0,0,0) => -4 // 1, (0,0,1) => +1 // 3, (0,1,1) => +1 // 6))
end
 
"Return the stencil for the 6th order left-hand side Mehrstellen operator"
function LeftMehrstellen6()
    return FullStencil(Dict(
        (0,0,0) => -64 // 15,
        (0,0,1) => +7 // 15,
        (0,1,1) => +1 // 10,
        (1,1,1) => +1 // 30,
    ))
end

"Return the 4th order 3D Mehrstellen stencil"
function Mehrstellen4()
    return MehrstellenStencil(
        Dict(
            (0,0,0) => -4 // 1, 
            (0,0,1) => +1 // 3, 
            (0,1,1) => +1 // 6
        ),
        Dict(
            (0,0,0) => +1 // 2,
            (0,0,1) => +1 // 12
        )
    )
end

"Return the 6th order 3D Mehrstellen stencil"
function Mehrstellen6()
    return MehrstellenStencil(
        Dict(
            (0,0,0) => -64 // 15,
            (0,0,1) => +7 // 15,
            (0,1,1) => +1 // 10,
            (1,1,1) => +1 // 30,
        ),
        Dict(
            (0,0,0) => +67 // 120,
            (0,0,1) => +1 // 18,
            (0,0,2) => -1 // 240,
            (0,1,1) => +1 // 90,
        )
    )
end

function StandardDifference(order::Int, ::Val{N}) where {N}
    @assert(order in (2, 4, 6, 8), "StandardDifference supports orders 2, 4, 6, and 8 (given $(order))")
    StencilType = SplitStencil{N, Rational{Int64}}
    if order == 2
        return StencilType([1//1])
    elseif order == 4
        return StencilType([4//3, -1//12])
    elseif order == 6
        return StencilType([3//2, -3//20, 1//90])
    else
        return StencilType([8//5, -1//5, 8//315, -1//560])
    end
end

"Return a 2D dimension-split stencil of given order and minimal width"
StandardDifference2D(order::Int) = StandardDifference(order, Val(2))

"Return a 3D dimension-split stencil of given order and minimal width"
StandardDifference3D(order::Int) = StandardDifference(order, Val(3))