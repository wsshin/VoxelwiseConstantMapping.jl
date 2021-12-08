export integral, integrate

const RTOL_DEFAULT = 1e-4

integral(args...; kwargs...) = integrate(args...; kwargs...).val

function integrate(vc::VoxelwiseConstant{K,AN}) where {K,AN}
    CI = CartesianIndices(size(vc))
    val = zero(eltype(AN))
    for ci = CI
        ∆lvxl_ci = map((∆lvxlₖ,iₖ) -> ∆lvxlₖ[iₖ], vc.∆lvxl, ci.I)
        val += vc.val[ci] * prod(∆lvxl_ci)
    end

    return (val=val, err=0.0)
end

function integrate(vc::VoxelwiseConstant,
                   dom;  # dom(x) returns true if x is inside domain; false otherwise
                   kwargs...)
    val, err = hcubature(x -> vc(x) * dom(x), bounds(vc)...; rtol=RTOL_DEFAULT, kwargs...)
    return (val=val, err=err)
end
