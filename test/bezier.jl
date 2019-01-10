using NURBS
using Base.Test

@testset "Bezier Curves" begin
    
    @testset "Bezier Curve" begin
        b = [1.0 2 4 3; 1 3 3 1; 0 0 0 0]
        ans1 = [1.0 1.56481 2.18519 2.75 3.14815 3.26852 3.0; 1.0 1.83333 2.33333  2.5 2.33333 1.83333 1.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0]
        ans2 = [[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7]]
        @test typeof(bezier(4, b, 7)[1]) == Array{Float64,2}
        @test typeof(bezier(4, b, 7)[2]) == Array{Array{Int64,1},1}
        @test isapprox(bezier(4, b, 7)[1], ans1, atol=1e-5)
        @test bezier(4, b, 7)[2] == ans2
    end

    @testset "Bezier Curves with derivatives" begin

    end
    
    @testset "Bernstein Basis" begin
        ans1 = [
            [1.0, 0.614, 0.275, 0.125, 0.042, 0.003, 0.0],
            [0.0, 0.325, 0.444, 0.375, 0.239, 0.058, 0.0],
            [0.0, 0.058, 0.239, 0.375, 0.444, 0.325, 0.0],
            [0.0, 0.003, 0.042, 0.125, 0.275, 0.614, 1.0],
        ]
        @test typeof(bern_basis(1, 1, 0.5)) == Float64
        @test isapprox(bern_basis(3, 2, 0.5), 0.375, atol = 1e-3)
        @test isapprox(map(y -> map(x -> bern_basis(3, y, x), [0, 0.15, 0.35, 0.5, 0.65, 0.85, 1]),[0, 1, 2, 3]), ans1, atol=1e-2)
    end
end
