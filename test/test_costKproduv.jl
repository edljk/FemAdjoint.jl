using Test, JLD2, UnicodePlots
using LinearAlgebra, ForwardDiff, Arpack
using SparseArrays, SparseDiffTools, Symbolics 

# load mesh
meshfile = "$(@__DIR__)/data/squaremesh.jld2"
JLD2.@load(meshfile, p, t, jacK, jacM, ps, ts, jacKs, jacMs)
#p, t, jacK, jacM =  ps, ts, jacKs, jacMs

np, nt = size(p, 1), size(t, 1)
ε = 1e-6
x, dp = p, 2 * rand(size(p)...) .- 1
xp, xm = p .+ ε * dp, p .- ε * dp
u, v = 2 * rand(np) .- 1, 2 * rand(np) .- 1

# finite differences
f = x -> FemAdjoint.costKproduv(x, t, u, v)
fp, fm = f(xp), f(xm)
g =  FemAdjoint.∇costKproduv(x, t, u, v, jacK = jacK, jacM = jacM)
println("")
println(fp)
println(fm)
println(" ")
println(dot(g, dp))
println((fp - fm) / (2 * ε))
println("")
@testset "test eval gradient dot function for costnormU" begin 
    @test abs(dot(g, dp) - (fp - fm) / (2 * ε)) < 1e-4
end