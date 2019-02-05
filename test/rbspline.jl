using NURBS
using Base.Test

@testset "RB-Spline curve generator" begin
   
    @testset "rbasis" begin
        x = knot(5,4)
        h = [1, 1, 0.25, 1, 1]
        @test typeof(rbasis(5,4,0.5,x,h)) == Array{Float64,1}
        @test isapprox(rbasis(5,4,0.5,x,h),[0.153846, 0.730769, 0.0769231, 0.0384615, 0.0],atol=1e-5)
    end
    
    @testset "rbspline" begin
        b = [0 1 5/2 4 5; 1 2 0 2 0; 0 0 0 0 0]
        h = [1, 1, 0.25, 1, 1]
        @test typeof(rbspline(5,3,10,b,h)[1]) == Array{Float64,2}
        @test typeof(rbspline(5,3,10,b,h)[2]) == Array{Array{Int64,1},1}
        @test isapprox(rbspline(5,3,10,b,h)[1],[0.0 0.557971 0.966667 1.3 1.95455 3.04545 3.7 4.03333 4.44203 5.0; 1.0 1.50725 1.73333 1.6 1.21212 1.21212 1.6 1.6 1.04348 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0],atol=1e-5)
    end
    
    @testset "rbsplinu" begin
        b = [0 1 5/2 4 5; 1 2 0 2 0; 0 0 0 0 0]
        h = [1, 1, 0.25, 1, 1]
        @test typeof(rbsplinu(5,3,10,b,h)[1]) == Array{Float64,2}
        @test typeof(rbsplinu(5,3,10,b,h)[2]) == Array{Array{Int64,1},1}
        @test isapprox(rbsplinu(5,3,10,b,h)[1],[0.5 0.789855 1.03333 1.3 1.95455 3.04545 3.7 3.96667 4.21014 4.75; 1.5 1.73913 1.8 1.6 1.21212 1.21212 1.6 1.73333 1.50725 0.5; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0],atol=1e-4)
    end
    
end