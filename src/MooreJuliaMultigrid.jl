module MooreJuliaMultigrid

using LinearAlgebra

"""
Run the Jacobi method on matrix system Ax=b for a total of numiters iterations.

    Start with an initial guess of x0
TODO: Change this to a weighted Jacobi - can use UniformScaling Matrix to simplify the weighted part
"""
function Jacobi(A,b,x0,numiters=1)
    Dinv = inv(Diagonal(diag(A)))
    R = Dinv*(A - Diagonal(diag(A)))
    xnew = copy(x0)
    for i = 1:numiters
        xnew = Dinv*b - R*xnew
    end
    return xnew
end

export Jacobi

end # module
