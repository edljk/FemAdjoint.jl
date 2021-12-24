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
