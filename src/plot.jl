#= to generate a square mesh 
nx, ny = 50, 50
X = [x for x ∈  range(-1, stop = 1., length = nx), y = 1:ny]
Y = [y for x = 1:nx, y ∈  range(-1, stop = 1., length = ny)]
p = hcat(X[:], Y[:])
t = delaunay(p)

nx, ny = 10, 10
Xs = [x for x ∈  range(-1, stop = 1., length = nx), y = 1:ny]
Ys = [y for x = 1:nx, y ∈  range(-1, stop = 1., length = ny)]
ps = hcat(Xs[:], Ys[:])
ts = delaunay(ps)

@save "/home/oudet/.julia/dev/FemAdjoint/test/data/squaremesh.jld" p t ps ts

Z = sin.(2 * X * pi) .* cos.(3 * Y * pi)
heatmap(Z, xfact=.1, yfact=.1, xoffset=-1.5, colormap = :inferno)
=#
