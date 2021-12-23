"""
   X, Y, p, t = gensquaremesh(nx::Int64,  ny::Int64 = 50)

to generate a square mesh 
"""
function gensquaremesh(nx::Int64,  ny::Int64 = 50)
    X = [x for x ∈  range(-1, stop = 1., length = nx), y = 1:ny]
    Y = [y for x = 1:nx, y ∈  range(-1, stop = 1., length = ny)]
    p = hcat(X[:], Y[:])
    t = delaunay(p)
    
    nx, ny = 10, 10
    Xs = [x for x ∈  range(-1, stop = 1., length = nx), y = 1:ny]
    Ys = [y for x = 1:nx, y ∈  range(-1, stop = 1., length = ny)]
    ps = hcat(Xs[:], Ys[:])
    ts = delaunay(ps)
    # save display
    @save "/home/oudet/.julia/dev/FemAdjoint/test/data/squaremesh.jld2" p t ps ts
    Z = sin.(2 * X * pi) .* cos.(3 * Y * pi)
    show(heatmap(Z, xfact=.1, yfact=.1, xoffset=-1.5, colormap = :inferno))
    return nothing
end

function genjac(meshfile::String = "diskmesh.jld2";
                addstr::String = "")
    meshfile = "/home/oudet/.julia/dev/FemAdjoint/test/data/" * meshfile 
    dictjld = FileIO.load(meshfile)
    p, t = dictjld["p" * addstr], dictjld["t" * addstr]
    if ("jacK" * addstr) ∉ keys(dictjld) # regen jac data
        y = similar(FemAdjoint.assembKM_P12D(p, t)[1])
        fs(y, x) = FemAdjoint.assembK_P12D_inplace(x, t, y)
        println("recompute sparsity pattern")
        @time sparsity_pattern = Symbolics.jacobian_sparsity(fs, y, p[:])
        jac = Float64.(sparse(sparsity_pattern))
        dictjld["jacK" * addstr] = jac
        FileIO.save(meshfile, dictjld)
        display(dictjld)
        println("end export sparsity pattern")
    end
    return nothing
end