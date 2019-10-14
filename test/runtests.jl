using MooreJuliaMultigrid
using LinearAlgebra
using Test

@testset "MooreJuliaMultigrid.jl" begin
    # Matrix system to test
    n = 128
    h = 1/(n-1)
    A = Tridiagonal(-1*ones(n-2),2*ones(n-1),-1*ones(n-2))/h^2
    b = zeros(n-1)
    x0 = [.2*sin(3*pi*x/n)+.7*sin(4*pi*x/n)+.4*sin(10*pi*x/n) for x = 1:(n-1)]
    prolongOp = BuildLinearInterpOp(n-1)
    restrictOp = 0.5*prolongOp'

    @test norm(Jacobi(A,b,x0,numiters=10,weight=2/3)) < 10^(-2)
    @test norm(BuildLinearInterpOp(7) - 0.5*[1 0 0;
                                         2 0 0;
                                         1 1 0;
                                         0 2 0;
                                         0 1 1;
                                         0 0 2;
                                         0 0 1]) == 0
    @test norm(Vcycle(A,b,x0,Jacobi,(numiters=2,weight=2/3),(numiters=2,weight=2/3),prolongOp,restrictOp)) < 10^(-2)

end
