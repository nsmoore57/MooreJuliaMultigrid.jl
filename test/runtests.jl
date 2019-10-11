using MooreJuliaMultigrid
using LinearAlgebra
using Test

@testset "MooreJuliaMultigrid.jl" begin
    # Write your own tests here.
    n = 10
    h = 1/n
    d = 2*ones(n)
    dl = -1*ones(n-1)
    A = Tridiagonal(dl,d,dl)/h^2
    b = zeros(n)
    x0 = [.2*sin(3*x/n)+.7*sin(4*x/n)+.4*sin(10*x/n) for x = 1:n]
    @test norm(Jacobi(A,b,x0,200,2/3)) < 0.001
end
