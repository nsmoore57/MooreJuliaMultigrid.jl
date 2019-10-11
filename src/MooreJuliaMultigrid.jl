module MooreJuliaMultigrid

using LinearAlgebra

"""
Run the Jacobi method on matrix system Ax=b for a total of numiters iterations.

    Start with an initial guess of x0
"""
function Jacobi(A,b,x0,numiters=1,weight=1)
    wDinv = weight*inv(Diagonal(diag(A)))
    R = UniformScaling(1-weight) + wDinv*(A - Diagonal(diag(A)))
    xnew = copy(x0)
    for i = 1:numiters
        xnew = wDinv*b - R*xnew
    end
    return xnew
end

export Jacobi

end # module
