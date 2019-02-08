using NURBS
using Base.Test

@testset "Rbsurfaces Building" begin
    
	vec = [-100.0 -100.0 -100.0 -100.0 -100.0 -50.0 -50.0 -50.0 -50.0 -50.0 0.0 0.0 0.0 0.0 0.0 50.0 50.0 50.0 50.0 50.0 100.0 100.0 100.0 100.0 100.0;-100.0 -50.0 0.0 50.0 100.0 -100.0 -50.0 0.0 50.0 100.0 -100.0 -50.0 0.0 50.0 100.0 -100.0 -50.0 0.0 50.0 100.0 -100.0 -50.0 0.0 50.0 100.0;0.0 0.0 0.0 0.0 0.0 0.0 25.0 50.0 25.0 0.0 0.0 25.0 50.0 25.0 0.0 0.0 25.0 150.0 25.0 0.0 0.0 0.0 0.0 0.0 0.0; 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0]

	res = [-100.0 -100.0 -100.0 -100.0 -25.9259 -25.9259 -25.9259 -25.9259 25.9259 25.9259 25.9259 25.9259 100.0 100.0 100.0 100.0; -100.0 -25.9259 25.9259 100.0 -100.0 -25.9259 25.9259 100.0 -100.0 -25.9259 25.9259 100.0 -100.0 -25.9259 25.9259 100.0; 0.0 0.0 0.0 0.0 0.0 34.8422 34.8422 0.0 0.0 51.3032 51.3032 0.0 0.0 0.0 0.0 0.0]

	@testset "rational B-spline surface using an open uniform knot vector" begin
		@test typeof(rbspsurf(vec,4,4,5,5,100,100)) == Tuple{Array{Float64, 2}, Array{Array{Int64, 1}, 1}, Array{Array{Int64, 1}, 1}}
		@test isapprox(rbspsurf(vec,4,4,5,5,4,4)[1], res; atol=1e-3)
	end
	
	@testset "matrix version" begin
		@test typeof(matrixnurbs(vec,4,4,5,5,100,100)) == Tuple{Array{Float64, 2}, Array{Array{Int64, 1}, 1}, Array{Array{Int64, 1}, 1}}
		@test typeof(matrixnurbs(vec,4,4,5,5,4,4)[1], res; atol=1e-3)
	end	
end