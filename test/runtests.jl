#using FemAdjoint
using Test, JLD2
using ForwardDiff

# load mesh and data for tests
meshfile = "$(@__DIR__)/data/diskmesh.jld"
JLD2.@load(meshfile, p, t)
np, nt = size(p, 1), size(t, 1)
ε = 1e-5
u, v = 2 * rand(np) .- 1, 2 * rand(np) .- 1
x, dp = p, 2 * rand(size(p)...) .- 1
xp, xm = p .+ ε * dp, p .- ε * dp

@testset "mesh data" begin   
    @test size(p, 1) == 560
    @test size(t, 1) == 1045 
end

@testset "first assembly" begin   
    SK, SM = FemAdjoint.assembKM_P12D(p, t)
    IK, JK = FemAdjoint.indKM_sparse(t)
    IM, JM = FemAdjoint.indKM_sparse(t)
    @test length(IK) == length(JK)
    @test length(IK) == length(SK)
    @test length(IM) == length(JM)
    @test length(IM) == length(SM)
end
K = sparse(IK, JK, SK, np, np)
M = sparse(IM, JM, SM, np, np)

@testset "test eval function" begin   
    u, v = 2 * rand(np) .- 1, 2 * rand(np) .- 1
    @test abs(dot(K * u, v) - FemAdjoint.costKproduv(p, t, u, v)) < 1e-8
end

f = x -> FemAdjoint.costKproduv(x, t, u, v)
fp, fm = f(xp), f(xm)
@time g = ForwardDiff.gradient(f, x)

@testset "test eval gradient dot function" begin 
    @test abs(dot(g, dp) - (fp - fm) / (2 * ε)) < 1e-4
end

display(K)
println("directional derivative $(abs(dot(g, dp)))")