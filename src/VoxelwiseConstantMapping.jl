module VoxelwiseConstantMapping

using StaticArrays
export VoxelwiseConstant

struct VoxelwiseConstant{K,AN<:AbstractArray{<:Number,K},VR<:AbstractVector{<:Real}}
    val::AN  # scalar values in voxels
    vxlbound::NTuple{K,VR}  # locations of faces of voxels including ghost locations; length(face[k]) = size(val,k) + 1

    function VoxelwiseConstant{K,AN,VR}(val, vxlbound) where {K,AN,VR}
        size(val) .+ 1 == length.(vxlbound) ||
            @error "length.(vxlbound) = $(length.(vxlbound)) should be greater by 1 in each dimension than size(val) = $(size(val))."

        for k = 1:K
            issorted(vxlbound[k]) || @error "vxlbound[$k] = $(vxlbound[k]) should be sorted."
        end

        return new(val, vxlbound)
     end
end

VoxelwiseConstant(val::AN, vxlbound::NTuple{K,VR}) where {K,AN<:AbstractArray{<:Number,K},VR<:AbstractVector{<:Real}} =
    VoxelwiseConstant{K,AN,VR}(val, vxlbound)

bounds(vc::VoxelwiseConstant) = (SVector(minimum.(vc.vxlbound)), SVector(maximum.(vc.vxlbound)))


(vc::VoxelwiseConstant{K})(pt::AbstractVector{<:Real}) where {K} = vc(SVector{K}(pt))
function (vc::VoxelwiseConstant{K})(pt::SVector{K,<:Real}) where {K}
    vxlbound = vc.vxlbound
    for k = 1:K
        vxlbound[k][1] ≤ pt[k] ≤ vxlbound[k][end] || @error "pt[$k] = $(pt[k]) should be between vxlbound[$k][1] = $(vxlbound[k][1]) and vxlbound[k][end] = $(vxlbound[k][end])."
    end

    # Find the pt-containing voxel with an inclusive negative-end boundary and exclusive
    # positive-end boundary.
    indₛ = ntuple(k -> findlast(w -> w≤pt[k], vxlbound[k][1:end-1]), Val(K))  # end-1 to ensure voxel's positive-end boundary is exclusive
    at_bound = ntuple(k -> (vxlbound[k][indₛ[k]]==pt[k] && indₛ[k]≥2), Val(K))
    indₑ = indₛ .+ at_bound
    CI = CartesianIndices(map((nᵢ,nₑ) -> nᵢ:nₑ, indₛ, indₑ))  # CartesianIndices{K}

    val = 0.0
    for ci = CI
        val += vc.val[ci]
    end
    val /= length(CI)

    return val
end

end # module
