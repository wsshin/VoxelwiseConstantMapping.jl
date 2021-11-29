using VoxelwiseConstantMapping
using Statistics: mean
using Test

@testset "VoxelwiseConstantMapping" begin
    N = (10,10)
    val = rand(N...)
    bnd_min, bnd_max = -5, 5
    vxlbound = (bnd_min:bnd_max, bnd_min:bnd_max)
    vc = VoxelwiseConstant(val, vxlbound)

    @test bounds(vc) == ([bnd_min,bnd_min], [bnd_max,bnd_max])
    @test vc([-0.5, -0.5]) == val[5, 5]
    @test vc([0, -0.5]) == mean(val[5:6, 5])
    @test vc([-0.5, 0]) == mean(val[5, 5:6])
    @test vc([0, 0]) == mean(val[5:6, 5:6])

end  # @testset "VoxelwiseConstantMapping"
