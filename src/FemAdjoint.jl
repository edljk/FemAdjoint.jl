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
    K = sparse(IK, JK, SK, np, np)
    M = sparse(IM, JM, SM, np, np)
    F = M * ones(np)
    U = K \ M
    # cost 
    val = sum(U[k] ^ 2 for k = 1:np)
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
    K = sparse(IK, JK, SK, np, np)
    M = sparse(IM, JM, SM, np, np)
    F = M * ones(np)
    U = K \ M
    # jacobian structures
    if length(jacK) == 0
        y = similar(FemAdjoint.assembKM_P12D(p, t)[1])
        fs(y, x) = FemAdjoint.assembK_P12D_inplace(x, t, y)
        sparsity_pattern = Symbolics.jacobian_sparsity(fs, y, p[:])
        jacK = Float64.(sparse(sparsity_pattern))
    end
    if length(jacM) == 0
        y = similar(FemAdjoint.assembKM_P12D(p, t)[2])
        fs(y, x) = FemAdjoint.assembM_P12D_inplace(x, t, y)
        sparsity_pattern = Symbolics.jacobian_sparsity(fs, y, p[:])
        jacM = Float64.(sparse(sparsity_pattern))
    end
    # adjoint
    JFU = 2 * U 
    H = - K \ JFU
    val = sum(U[k] ^ 2 for k = 1:np)
    return val
end

end # module