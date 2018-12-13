using NURBS
using Base.Test

@testset "B-spline Building" begin
        
    @testset "B-spline curve using matrix methods and a periodic uniform knot vector" begin

        vector = [2.0 4.0 4.0 4.0 2.0 0.0 0.0 0.0 2.0; 0.0 0.0 2.0 4.0 4.0 4.0 2.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]
	
        res = [4.0 3.66667 2.0 0.333333 0.0 0.333333 1.66667 2.33333 3.66667; 2.0 3.66667 4.0 3.66667 2.0 0.333333 0.0 0.0 0.333333; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]

        @test typeof(matpbspl(4, 9, 200,vector)) == Array{Float64,2}
        @test isapprox(matpbspl(4, 9, 1,vector), res; atol=1e-3)
        end
    
    @testset "B-spline periodic basis matrix" begin
        @test typeof(nmatrix(4)) == Tuple{Array{Float64,2}, Float64}
        @test nmatrix(4) == ([-1.0 3.0 -3.0 1.0; 3.0 -6.0 3.0 0.0; -3.0 0.0 3.0 0.0; 1.0 4.0 1.0 0.0], 0.16666666666666666)
        @test nmatrix(5) == ([1.0 -4.0 6.0 -4.0 1.0; -4.0 12.0 -12.0 4.0 0.0; 6.0 -6.0 -6.0 6.0 0.0; -4.0 -12.0 12.0 4.0 0.0; 1.0 11.0 11.0 1.0 0.0], 0.041666666666666664)
    end
end 
