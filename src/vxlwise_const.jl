export VoxelwiseConstant
export bounds

struct VoxelwiseConstant{K,AN<:AbsArrNumber{K}}
    val::AN  # scalar values in voxels
    lvxl::NTuple{K,VecFloat}  # locations of faces (not centers) of voxels
    ∆lvxl::NTuple{K,VecFloat}  # voxel edges


    function VoxelwiseConstant{K,AN}(val, lvxl) where {K,AN}
        size(val) .+ 1 == length.(lvxl) ||
            @error "length.(lvxl) = $(length.(lvxl)) should be greater by 1 in each dimension than size(val) = $(size(val))."

        for k = 1:K
            issorted(lvxl[k]) || @error "lvxl[$k] = $(lvxl[k]) should be sorted."
        end

        ∆lvxl = diff.(lvxl)
        @assert length.(∆lvxl)==size(val)

        return new(val, lvxl, ∆lvxl)
     end
end

VoxelwiseConstant(val::AN, lvxl::NTuple{K,AbsVecReal}) where {K,AN<:AbsArrNumber{K}} =
    VoxelwiseConstant{K,AN}(val, lvxl)

Base.size(vc::VoxelwiseConstant) = size(vc.val)
bounds(vc::VoxelwiseConstant) = (SVec(minimum.(vc.lvxl)), SVec(maximum.(vc.lvxl)))

(vc::VoxelwiseConstant{K})(x::AbsVecReal) where {K} = vc(SFloat{K}(x))

function (vc::VoxelwiseConstant{K,AN})(x::SReal{K}) where {K,AN<:AbsArrNumber{K}}
    lvxl = vc.lvxl
    for k = 1:K
        lvxl[k][1] ≤ x[k] ≤ lvxl[k][end] || @error "x[$k] = $(x[k]) should be between lvxl[$k][1] = $(lvxl[k][1]) and lvxl[k][end] = $(lvxl[k][end])."
    end

    # Find the x-containing voxel with an inclusive negative-end boundary and exclusive
    # positive-end boundary.
    ind_rng = MVec{K,UnitRange{Int64}}(undef)
    for k = 1:K
        xₖ = x[k]
        lvxlₖ = lvxl[k]

        indₑ = findlast(w -> w≤xₖ, @view(lvxlₖ[1:end-1]))

        at_bound = lvxlₖ[indₑ]==xₖ && indₑ≥2
        indₛ = indₑ - at_bound

        ind_rng[k] = indₛ:indₑ
    end
    CI = CartesianIndices(ind_rng.data)  # CartesianIndices{K}

    return mean(@view(vc.val[CI]))
end
