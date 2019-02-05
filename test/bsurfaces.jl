using NURBS
using Base.Test

@testset "bsurface Building" begin
    
	arr = [-100.0 -100.0 -100.0 -100.0 -100.0 -50.0 -50.0 -50.0 -50.0 -50.0 0.0 0.0 0.0 0.0 0.0 50.0 50.0 50.0 50.0 50.0 100.0 100.0 100.0 100.0 100.0;-100.0 -50.0 0.0 50.0 100.0 -100.0 -50.0 0.0 50.0 100.0 -100.0 -50.0 0.0 50.0 100.0 -100.0 -50.0 0.0 50.0 100.0 -100.0 -50.0 0.0 50.0 100.0;0.0 0.0 0.0 0.0 0.0 0.0 25.0 50.0 25.0 0.0 0.0 25.0 50.0 25.0 0.0 0.0 25.0 150.0 25.0 0.0 0.0 0.0 0.0 0.0 0.0]

	res = [-100.0 -100.0 -100.0 -100.0 -25.9259 -25.9259 -25.9259 -25.9259 25.9259 25.9259 25.9259 25.9259 100.0 100.0 100.0 100.0; -100.0 -25.9259 25.9259 100.0 -100.0 -25.9259 25.9259 100.0 -100.0 -25.9259 25.9259 100.0 -100.0 -25.9259 25.9259 100.0; 0.0 0.0 0.0 0.0 0.0 34.8422 34.8422 0.0 0.0 51.3032 51.3032 0.0 0.0 0.0 0.0 0.0]

    resu = [
        -50.0 -50.0 -50.0 -50.0 -16.6667 -16.6667 -16.6667  -16.6667 16.6667 16.6667 16.6667 16.6667 50.0 50.0 50.0 50.0; 
        -50.0 -16.6667 16.6667 50.0 -50.0 -16.6667 16.6667 50.0 -50.0 -16.6667 16.6667 50.0 -50.0 -16.6667 16.6667 50.0;
        20.8333 32.6646 32.6646 20.8333 25.6687 41.7905 41.7905 25.6687 31.0185 60.2176 60.2176 31.0185 31.9444 70.9362 70.9362 31.9444
    ]

    @testset "B-spline surface using open uniform knot vectors" begin
        @test typeof(bsplsurf(arr, 4, 4, 5, 5, 100, 100)) == Tuple{Array{Float64, 2}, Array{Array{Int64, 1}, 1}, Array{Array{Int64, 1}, 1}}
        @test isapprox(bsplsurf(arr, 4, 4, 5, 5, 4, 4)[1], res; atol = 1e-3)
    end

    @testset "B-spline surface using periodic uniform knot vectors" begin
        @test typeof(bsplsurfu(arr, 4, 4, 5, 5, 100, 100)) == Tuple{Array{Float64, 2}, Array{Array{Int64, 1}, 1}, Array{Array{Int64, 1}, 1}}
        @test isapprox(bsplsurfu(arr, 4, 4, 5, 5, 4, 4)[1], resu; atol = 1e-3)
    end
        
    @testset "B-spline surface and its derivatives using open uniform knot vectors" begin
        @test typeof(dbsurf(arr, 4, 4, 5, 5, 100, 100)) == Tuple{Array{Float64,2}, Array{Float64,2}, Array{Float64,2}, Array{Float64,2}, Array{Float64,2}, Array{Float64,2}, Array{Array{Int64, 1}, 1}, Array{Array{Int64, 1}, 1}}
        @test isapprox(dbsurf(arr, 4, 4, 5, 5, 4, 4)[1], res; atol = 1e-3)
    end

    @testset "B-Spline surface over open knot comparison" begin
        @test bsplsurf(arr, 4, 4, 5, 5, 100, 100)[1] == dbsurf(arr, 4, 4, 5, 5, 100, 100)[1]
        @test bsplsurf(arr, 4, 4, 5, 5, 100, 100)[2] == dbsurf(arr, 4, 4, 5, 5, 100, 100)[7]
        @test bsplsurf(arr, 4, 4, 5, 5, 100, 100)[3] == dbsurf(arr, 4, 4, 5, 5, 100, 100)[8]
    end
end
