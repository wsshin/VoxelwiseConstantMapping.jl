@testset "integral, constant" begin
    N = (10,20)
    v₀ = rand()
    val = fill(v₀, N)

    p₀ = [rand(), rand()]
    r = 5.0
    xₙ, xₚ = p₀[1]-r, p₀[1]+r
    yₙ, yₚ = p₀[2]-r, p₀[2]+r

    xbound = range(xₙ, xₚ, length=N[1]+1)
    ybound = range(yₙ, yₚ, length=N[2]+1)
    vxlbound = (xbound, ybound)

    vc = VoxelwiseConstant(val, vxlbound)

    @test integral(vc) ≈ v₀ * (xₚ-xₙ) * (yₚ-yₙ)
    @test (s = integrate(vc, p->norm(p .- p₀)≤r); abs(s.val - v₀ * π * r^2) ≤ s.err)
end  # @testset "integral, constant"

@testset "integral, non-constant" begin
    N = (10,20)
    val = rand(N...)

    p₀ = [rand(), rand()]
    r = 5.0
    xₙ, xₚ = p₀[1]-r, p₀[1]+r
    yₙ, yₚ = p₀[2]-r, p₀[2]+r

    xbound = range(xₙ, xₚ, length=N[1]+1)
    ybound = range(yₙ, yₚ, length=N[2]+1)
    vxlbound = (xbound, ybound)

    ∆x = (xₚ-xₙ) / N[1]
    ∆y = (yₚ-yₙ) / N[2]

    vc = VoxelwiseConstant(val, vxlbound)

    @test integral(vc) ≈ sum(val) * ∆x * ∆y
end  # @testset "integral, non-constant"
