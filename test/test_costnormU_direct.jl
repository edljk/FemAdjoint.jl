using Test, JLD2, UnicodePlots
using LinearAlgebra, ForwardDiff, Arpack
using SparseArrays, SparseDiffTools, Symbolics 

# load mesh
meshfile = "$(@__DIR__)/data/squaremesh.jld2"
JLD2.@load(meshfile, p, t, ps, ts)
#p, t =  ps, ts

np, nt = size(p, 1), size(t, 1)
ε = 1e-4
x, dp = p, 2 * rand(size(p)...) .- 1
Ib = unique(FemAdjoint.btri(t)[:])
dp[Ib, :] .= 0.# fix Dirichlet points
xp, xm = p .+ ε * dp, p .- ε * dp

# finite differences
f = x -> FemAdjoint.costnormU(x, t)
fp, fm = f(xp), f(xm)
@time g = FemAdjoint.∇costnormU_direct(x, t, plotsolutions = true)
@time g = FemAdjoint.∇costnormU_direct(x, t, plotsolutions = false)
println("")
println(fp)
println(fm)
println(" ")
println((fp - fm) / (2 * ε))
println(dot(g, dp))
println("")
@testset "test eval gradient dot function for costnormU" begin 
    @test abs(dot(g, dp) - (fp - fm) / (2 * ε)) < 1e-4
end