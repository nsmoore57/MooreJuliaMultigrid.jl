module MooreJuliaMultigrid

using LinearAlgebra
using SparseArrays

"""
Run the Weighted Jacobi method on matrix system Ax=b for a total of numiters
    iterations, starting with an initial guess of x0
"""
function Jacobi(A,b,x0; numiters=1,weight=1)
    wDinv = weight*inv(Diagonal(diag(A)))
    R = UniformScaling(1-weight) + wDinv*(A - Diagonal(diag(A)))
    for i = 1:numiters
        x0 = wDinv*b - R*x0
    end
    return x0
end

"""
V-cycle method to solve the matrix system Ax=b
Inputs:
    A, b, x0 : A matrix, RHS vector, initial guess
    relaxFunc : Function to complete the relaxation
      -- argument form : relaxFunc(A,b,x0; named options)
    preRelaxOpts : Options for relaxFunc before coarse-grid transfer
      -- named tuple with named options for relaxFunc
    postRelaxOpts : Options for relaxFunc after coarse-grid correction
      -- named tuple with named options for relaxFunc
    prolongOp : Prolongation (interpolation) operator
      -- matrix expected - may change
    restrictOp : Restriction operator
      -- matrix expected - may change
"""
function Vcycle(A,b,x0,relaxFunc,preRelaxOpts, postRelaxOpts, prolongOp, restrictOp)
    # First relax numprerelax times
    if preRelaxOpts[:numiters] > 0
        x0 = relaxFunc(A,b,x0;preRelaxOpts...)
    end

    # Transfer residual to coarser grid
    r = restrictOp*(b - A*x0)

    # Transfer A to coarser grid
    Acoarse = restrictOp*A*prolongOp

    # Solve the coarse problem
    ecourse = Acoarse \ r

    # Interp ecourse to fine grid
    e = prolongOp*ecourse

    # Correct approximation with the residual
    x0 = x0 + e

    # Relax numpostrelax times
    if postRelaxOpts[:numiters] > 0
        x0 = relaxFunc(A,b,x0;postRelaxOpts...)
    end
    return x0
end

"""
BuildLinearInterpOp builds the linear interpoltion operator for transferring
    a vector of length fineSize to a vector of length (fineSize+1)/2 - 1.

    Currently expects only an odd length vector

    TODO: Add support for even length vectors
"""
function BuildLinearInterpOp(fineSize)
    numCols::Int = (fineSize+1)/2 - 1
    I = Vector{Int}(undef,3*numCols)
    J = Vector{Int}(undef,3*numCols)
    V = Vector{Float64}(undef,3*numCols)
    for i = 1:numCols
        I[3i-2:3i] = [2i-1; 2i; 2i+1]
        J[3i-2:3i] = i*ones(3,1)
        V[3i-2:3i] = [1; 2; 1]
    end
    return 0.5*sparse(I,J,V)
end

export Jacobi, Vcycle, BuildLinearInterpOp

end # module
