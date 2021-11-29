module VoxelwiseConstantMapping

using AbbreviatedTypes

export VoxelwiseConstant
export bounds

struct VoxelwiseConstant{K,AN<:AbsArrNumber{K}}
    val::AN  # scalar values in voxels
    vxlbound::NTuple{K,VecFloat}  # locations of faces of voxels

    function VoxelwiseConstant{K,AN}(val, vxlbound) where {K,AN}
        size(val) .+ 1 == length.(vxlbound) ||
            @error "length.(vxlbound) = $(length.(vxlbound)) should be greater by 1 in each dimension than size(val) = $(size(val))."

        for k = 1:K
            issorted(vxlbound[k]) || @error "vxlbound[$k] = $(vxlbound[k]) should be sorted."
        end

        return new(val, vxlbound)
     end
end

VoxelwiseConstant(val::AN, vxlbound::NTuple{K,AbsVecReal}) where {K,AN<:AbsArrNumber{K}} =
    VoxelwiseConstant{K,AN}(val, vxlbound)

bounds(vc::VoxelwiseConstant) = (SVec(minimum.(vc.vxlbound)), SVec(maximum.(vc.vxlbound)))

(vc::VoxelwiseConstant{K})(pt::AbsVecReal) where {K} = vc(SFloat{K}(pt))

function (vc::VoxelwiseConstant{K,AN})(pt::SReal{K}) where {K,AN<:AbsArrNumber{K}}
    vxlbound = vc.vxlbound
    for k = 1:K
        vxlbound[k][1] ≤ pt[k] ≤ vxlbound[k][end] || @error "pt[$k] = $(pt[k]) should be between vxlbound[$k][1] = $(vxlbound[k][1]) and vxlbound[k][end] = $(vxlbound[k][end])."
    end

    # Find the pt-containing voxel with an inclusive negative-end boundary and exclusive
    # positive-end boundary.
    indₑ = ntuple(k -> findlast(w -> w≤pt[k], vxlbound[k][1:end-1]), Val(K))  # end-1 to ensure voxel's positive-end boundary is exclusive
    at_bound = ntuple(k -> (vxlbound[k][indₑ[k]]==pt[k] && indₑ[k]≥2), Val(K))
    indₛ = indₑ .- at_bound
    CI = CartesianIndices(map((nₛ,nₑ) -> nₛ:nₑ, indₛ, indₑ))  # CartesianIndices{K}

    val = 0.0
    for ci = CI
        val += vc.val[ci]
    end
    val /= length(CI)

    return val
end

end # module
