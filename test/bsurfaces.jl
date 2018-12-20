using NURBS
using Base.Test

@testset "bsurface Building" begin
    
	arr = [-100.0 -100.0 -100.0 -100.0 -100.0 -50.0 -50.0 -50.0 -50.0 -50.0 0.0 0.0 0.0 0.0 0.0 50.0 50.0 50.0 50.0 50.0 100.0 100.0 100.0 100.0 100.0;-100.0 -50.0 0.0 50.0 100.0 -100.0 -50.0 0.0 50.0 100.0 -100.0 -50.0 0.0 50.0 100.0 -100.0 -50.0 0.0 50.0 100.0 -100.0 -50.0 0.0 50.0 100.0;0.0 0.0 0.0 0.0 0.0 0.0 25.0 50.0 25.0 0.0 0.0 25.0 50.0 25.0 0.0 0.0 25.0 150.0 25.0 0.0 0.0 0.0 0.0 0.0 0.0]

	res = [-100.0 -100.0 -100.0 -100.0 -25.9259 -25.9259 -25.9259 -25.9259 25.9259 25.9259 25.9259 25.9259 100.0 100.0 100.0 100.0; -100.0 -25.9259 25.9259 100.0 -100.0 -25.9259 25.9259 100.0 -100.0 -25.9259 25.9259 100.0 -100.0 -25.9259 25.9259 100.0; 0.0 0.0 0.0 0.0 0.0 34.8422 34.8422 0.0 0.0 51.3032 51.3032 0.0 0.0 0.0 0.0 0.0]

    @testset "B-spline surface using open uniform knot vectors" begin
        @test typeof(bsplsurf(arr,4,4,5,5,100,100)) == Array{Float64,2}
        @test isapprox(bsplsurf(arr,4,4,5,5,4,4), res; atol=1e-3)
    end
        
    @testset "B-spline surface and its derivatives using open uniform knot vectors" begin
        @test typeof(dbsurf(arr,4,4,5,5,100,100)) == Tuple{Array{Float64,2}, Array{Float64,2}, Array{Float64,2}, Array{Float64,2}, Array{Float64,2}, Array{Float64,2}}
        @test isapprox(dbsurf(arr,4,4,5,5,4,4)[1], res; atol=1e-3)
    end
end
