module FemAdjoint

using LinearAlgebra, StatsBase, GroupSlices
using ForwardDiff, SparseDiffTools, SparseArrays, Symbolics
using FileIO, JLD2, UnicodePlots

include("assembly.jl")
include("plot.jl")

function costKproduv(pin, t, u, v)
    p = if typeof(pin) <: Vector
        reshape(pin, length(pin) ÷ 2, 2)
    elseif typeof(pin) <: Matrix
        pin
    end
    np = size(p, 1)
    I, J = indKM_sparse(t)
    SK, SM = FemAdjoint.assembKM_P12D(p, t)
    val = 0.
    for (i, j, s) ∈ zip(I, J, SK)
        val += u[i] * s * v[j]
    end
    return val
end
function ∇costKproduv(pin, t, u, v;
                      jacK::SparseMatrixCSC{Float64, Int64} = spzeros(0, 0),
                      jacM::SparseMatrixCSC{Float64, Int64} = spzeros(0, 0))
    p = if typeof(pin) <: Vector
        reshape(pin, length(pin) ÷ 2, 2)
    elseif typeof(pin) <: Matrix
        pin
    end
    np = size(p, 1)
    I, J = indKM_sparse(t)
    # jacobian structures
    fK(y, x) = FemAdjoint.assembK_P12D_inplace!(x, t, y)
    if length(jacK) == 0
        y = similar(FemAdjoint.assembKM_P12D(p, t)[1])
        sparsity_pattern = Symbolics.jacobian_sparsity(fK, y, p[:])
        jacK = Float64.(sparse(sparsity_pattern))
    end
    # evaluate jacobians 
    colorsK = matrix_colors(jacK)
    JK = forwarddiff_color_jacobian!(jacK, fK, p[:], colorvec = colorsK)
    # compute full gradient 
    g = zeros(np, 2)
    for (i, j, s) ∈ zip(findnz(JK)...)
        g[j] += u[J[i]] * s * v[I[i]]
    end
    return g
end
function costMproduv(pin, t, u, v)
    p = if typeof(pin) <: Vector
        reshape(pin, length(pin) ÷ 2, 2)
    elseif typeof(pin) <: Matrix
        pin
    end
    np = size(p, 1)
    I, J = indKM_sparse(t)
    SK, SM = FemAdjoint.assembKM_P12D(p, t)
    val = 0.
    for (i, j, s) ∈ zip(I, J, SM)
        val += u[i] * s * v[j]
    end
    return val
end
function ∇costMproduv(pin, t, u, v;
                      jacK::SparseMatrixCSC{Float64, Int64} = spzeros(0, 0),
                      jacM::SparseMatrixCSC{Float64, Int64} = spzeros(0, 0))
    p = if typeof(pin) <: Vector
        reshape(pin, length(pin) ÷ 2, 2)
    elseif typeof(pin) <: Matrix
        pin
    end
    np = size(p, 1)
    I, J = indKM_sparse(t)
    # jacobian structures
    fM(y, x) = FemAdjoint.assembM_P12D_inplace!(x, t, y)
    if length(jacM) == 0
        y = similar(FemAdjoint.assembKM_P12D(p, t)[2])
        sparsity_pattern = Symbolics.jacobian_sparsity(fK, y, p[:])
        jacK = Float64.(sparse(sparsity_pattern))
    end
    # evaluate jacobians 
    colorsM = matrix_colors(jacM)
    JM = forwarddiff_color_jacobian!(jacM, fM, p[:], colorvec = colorsM)
    # compute full gradient 
    g = zeros(np, 2)
    for (i, j, s) ∈ zip(findnz(JM)...)
        g[j] += u[J[i]] * s * v[I[i]]
    end
    return g
end
#-------------------------------------------------------------------------------
function costnormU(pin, t)
    p = if typeof(pin) <: Vector
        reshape(pin, length(pin) ÷ 2, 2)
    elseif typeof(pin) <: Matrix
        pin
    end
    np = size(p, 1)
    # assemble matrices
    I, J = indKM_sparse(t)
    SK, SM = FemAdjoint.assembKM_P12D(p, t)
    K = sparse(I, J, SK, np, np)
    Ib = unique(FemAdjoint.btri(t)[:])
    K += sparse(Ib, Ib, fill(1e8, length(Ib)), np, np) # Dirichlet conditions
    M = sparse(I, J, SM, np, np)
    F = M * ones(np)
    U = K \ F
    # cost 
    val = sum(U[k] ^ 2 for k = 1:np) / np
    return val
end

function ∇costnormU(pin, t;
                    jacK::SparseMatrixCSC{Float64, Int64} = spzeros(0, 0),
                    jacM::SparseMatrixCSC{Float64, Int64} = spzeros(0, 0))
    p = if typeof(pin) <: Vector
        reshape(pin, length(pin) ÷ 2, 2)
    elseif typeof(pin) <: Matrix
        pin
    end
    np = size(p, 1)
    # assemble matrices
    I, J = indKM_sparse(t)
    SK, SM = FemAdjoint.assembKM_P12D(p, t)
    K = sparse(I, J, SK, np, np)
    Ib = unique(FemAdjoint.btri(t)[:])
    K += sparse(Ib, Ib, fill(1e8, length(Ib)), np, np) # Dirichlet conditions
    M = sparse(I, J, SM, np, np)
    # compute state function
    F = M * ones(np)
    U = K \ F
    # plot
    plotunicode(U)
    # jacobian structures
    fK(y, x) = FemAdjoint.assembK_P12D_inplace!(x, t, y)
    fM(y, x) = FemAdjoint.assembM_P12D_inplace!(x, t, y)
    if length(jacK) == 0
        y = similar(FemAdjoint.assembKM_P12D(p, t)[1])
        sparsity_pattern = Symbolics.jacobian_sparsity(fK, y, p[:])
        jacK = Float64.(sparse(sparsity_pattern))
    end
    if length(jacM) == 0
        y = similar(FemAdjoint.assembKM_P12D(p, t)[2])
        sparsity_pattern = Symbolics.jacobian_sparsity(fM, y, p[:])
        jacM = Float64.(sparse(sparsity_pattern))
    end
    # evaluate jacobians 
    colorsK = matrix_colors(jacK)
    JK = forwarddiff_color_jacobian!(jacK, fK, p[:], colorvec = colorsK)
    colorsM = matrix_colors(jacM)
    JM = forwarddiff_color_jacobian!(jacM, fM, p[:], colorvec = colorsM)
    # adjoint function
    JFU = - 2 * U / np
    λ = K' \ JFU
    # compute full gradient 
    g = zeros(np, 2)
    for (i, j, s) ∈ zip(findnz(JK)...)
        if I[i] != J[i] && I[i] ∉ Ib
            g[j] += λ[J[i]] * s * U[I[i]]
        end
    end
    for (i, j, s) ∈ zip(findnz(JM)...)
        g[j] -= λ[J[i]] * s 
    end
    return g
end

end # module