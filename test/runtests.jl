#using FemAdjoint
using Test, JLD2, UnicodePlots
using LinearAlgebra, ForwardDiff, Arpack
using SparseArrays, SparseDiffTools, Symbolics 
using FemAdjoint

# load mesh and data for tests
meshfile = "$(@__DIR__)/data/diskmesh.jld2"
JLD2.@load(meshfile, p, t)
np, nt = size(p, 1), size(t, 1)
ε = 1e-5
u, v = 2 * rand(np) .- 1, 2 * rand(np) .- 1
x, dp = p, 2 * rand(size(p)...) .- 1
xp, xm = p .+ ε * dp, p .- ε * dp
SK, SM = FemAdjoint.assembKM_P12D(p, t)
IK, JK = FemAdjoint.indKM_sparse(t)
IM, JM = FemAdjoint.indKM_sparse(t)
K = sparse(IK, JK, SK, np, np)
M = sparse(IM, JM, SM, np, np)

@testset "mesh data" begin   
    @test size(p, 1) == 560
    @test size(t, 1) == 1045 
end

@testset "first assembly" begin   
    @test length(IK) == length(JK)
    @test length(IK) == length(SK)
    @test length(IM) == length(JM)
    @test length(IM) == length(SM)
end


@testset "test eval function" begin   
    u, v = 2 * rand(np) .- 1, 2 * rand(np) .- 1
    @test abs(dot(K * u, v) - FemAdjoint.costKproduv(p, t, u, v)) < 1e-8
end

f = x -> FemAdjoint.costKproduv(x, t, u, v)
fp, fm = f(xp), f(xm)
@time g = ForwardDiff.gradient(f, x)
@time g = ForwardDiff.gradient(f, x)

@testset "test eval gradient dot function" begin 
    @test abs(dot(g, dp) - (fp - fm) / (2 * ε)) < 1e-4
end
display(K)
print("\n")
println("directional derivative $(abs(dot(g, dp)))")

#-------------------------------------------------------------------------------
# square mesh / eigenvalue test
meshfile = "$(@__DIR__)/data/squaremesh.jld2"
JLD2.@load(meshfile, p, t, jacK, jacM, ps, ts, jacKs, jacMs)
np, nt = size(p, 1), size(t, 1)
SK, SM = FemAdjoint.assembKM_P12D(p, t)
IK, JK = FemAdjoint.indKM_sparse(t)
IM, JM = FemAdjoint.indKM_sparse(t)
K = sparse(IK, JK, SK, np, np)
M = sparse(IM, JM, SM, np, np)
Ib = unique(FemAdjoint.btri(t)[:])
K += sparse(Ib, Ib, fill(1e8, length(Ib)), np, np) # Dirichlet conditions

λ, U = eigs(K, M, nev = 10, which = :SM, tol = 0., maxiter = 30_000)
@testset "test eigenvalue" begin
    @test abs(λ[1] - pi ^ 2 / 2) < 1e-2
end
# plot
FemAdjoint.plotunicode(U[:, 4])

#-------------------------------------------------------------------------------
#p, t, jac = ps, ts, jacs
y = similar(FemAdjoint.assembKM_P12D(p, t)[1])
SK, SM = FemAdjoint.assembKM_P12D(p, t)
SKi, SMi = similar(SK), similar(SM)
FemAdjoint.assembK_P12D_inplace!(p, t, SKi)
FemAdjoint.assembM_P12D_inplace!(p, t, SMi)
@testset "test in place assemble" begin
    @test maximum(abs.(SK .- SKi)) < 1e-8
    @test maximum(abs.(SM .- SMi)) < 1e-8
end
fK(y, x) = FemAdjoint.assembK_P12D_inplace!(x, t, y)
# to speedup test use a precomputed jac matrix
#@time sparsity_pattern = Symbolics.jacobian_sparsity(fs, y, p[:])
#jacK = Float64.(sparse(sparsity_pattern))
colors = matrix_colors(jacK)
@time J = forwarddiff_color_jacobian!(jacK, fK, p[:], colorvec = colors)
@time J = forwarddiff_color_jacobian!(jacK, fK, p[:], colorvec = colors)
display(J)

#-------------------------------------------------------------------------------
# test adjoint approach for costKproduv cost function
include("test_costKproduv.jl")
#-------------------------------------------------------------------------------
# test adjoint approach for costnormU cost function
include("test_costnormU.jl")
#-------------------------------------------------------------------------------
# test direct adjoint approach for costnormU cost function (definitely simpler!!)
include("test_costnormU_direct.jl")