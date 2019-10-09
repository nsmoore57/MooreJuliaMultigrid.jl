using MooreJuliaMultigrid
using LinearAlgebra
using Test

@testset "MooreJuliaMultigrid.jl" begin
    # Write your own tests here.
    n = 100
    h = 1/n
    d = 2*ones(n)
    dl = -1*ones(n-1)
    A = Tridiagonal(dl,d,dl)/h^2
    b = zeros(n)
    x0 = rand(n)
    @test norm(Jacobi(A,b,x0,100000)) < 0.001
end
