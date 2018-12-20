using NURBS
using Base.Test

@testset "Knot Building" begin
    
    @testset "Open Uniform Knots" begin
        @test typeof(knot(5, 2)) == Array{Float64,1}
        @test knot(5, 2) == [0.0, 0.0, 1.0, 2.0, 3.0, 4.0, 4.0]
        @test knot(6, 3, false, 0.25) == [0.0, 0.0, 0.0, 0.25, 0.5, 0.75, 1.0, 1.0, 1.0]
        @test knot(5, 2, true) == [-2.0, -2.0, -1.0, 0.0, 1.0, 2.0, 2.0]
    end

    @testset "Open NonUniform Proportional Knots" begin
        b = [0. 2 4 6 8; 0 6 3 6 6; 0 0 0 0 0]
        @test typeof(knotc(5, 3, b)) == Array{Float64,1}
        @test isapprox(knotc(5, 3, b), [0.0, 0.0, 0.0, 1.45338, 2.38171, 3.0, 3.0, 3.0]; atol = 1e-6)
    end
    
    @testset "Periodic Uniform Knots" begin
        @test typeof(knotu(5, 2)) == Array{Float64,1}
        @test knotu(5,2) == [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
        @test knotu(5, 4, 0.25) == [0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0]
    end
end
