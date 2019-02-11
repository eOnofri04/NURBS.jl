using NURBS
using Base.Test

@testset "B-Spline curve generator" begin
   
    @testset "bspline" begin
        b = [1. 2 4 3; 1 3 3 1; 0 0 0 0]
        @test typeof(bspline(4,4,5,b)[1]) == Array{Float64,2}
        @test typeof(bspline(4,4,5,b)[2]) == Array{Array{Int64,1},1}
        @test isapprox(bspline(4,4,5,b)[1],[1.0 1.875 2.75 3.25 3.0; 1.0 2.125 2.5 2.125 1.0; 0.0 0.0 0.0 0.0 0.0],atol=1e-5)
    end
    
    @testset "bsplineu" begin
        b = [0. 3 6 9; 0 10 3 6; 0 0 0 0]
        @test typeof(bsplineu(4,3,5,b)[1]) == Array{Float64,2}
        @test typeof(bsplineu(4,3,5,b)[2]) == Array{Array{Int64,1},1}
        @test isapprox(bsplineu(4,3,5,b)[1],[1.5 3.0 4.5 6.0 7.5; 5.0 7.875 6.5 4.25 4.5; 0.0 0.0 0.0 0.0 0.0],atol=1e-5)
    end
    
    @testset "dbspline" begin
        b = [1. 2 4 3; 1 3 3 1; 0 0 0 0]
        @test typeof(dbspline(4,4,5,b)[1]) == Array{Float64,2}
        @test typeof(dbspline(4,4,5,b)[2]) == Array{Float64,2}
        @test typeof(dbspline(4,4,5,b)[3]) == Array{Float64,2}
        @test typeof(dbspline(4,4,5,b)[4]) == Array{Array{Int64,1},1}
        @test isapprox(dbspline(4,4,5,b)[1],[1.0 1.875 2.75 3.25 3.0; 1.0 2.125 2.5 2.125 1.0; 0.0 0.0 0.0 0.0 0.0],atol=1e-5)
        @test isapprox(dbspline(4,4,5,b)[2],[3.0 3.75 3.0 0.75 -3.0; 6.0 3.0 0.0 -3.0 -6.0; 0.0 0.0 0.0 0.0 0.0],atol=1e-5)
        @test isapprox(dbspline(4,4,5,b)[3],[6.0 0.0 -6.0 -12.0 -18.0; -12.0 -12.0 -12.0 -12.0 -12.0; 0.0 0.0 0.0 0.0 0.0],atol=1e-5)
    end
    
    @testset "dbsplineu" begin
        b = [0. 3 6 9; 0 10 3 6; 0 0 0 0]
        @test typeof(dbsplineu(4,3,5,b)[1]) == Array{Float64,2}
        @test typeof(dbsplineu(4,3,5,b)[2]) == Array{Float64,2}
        @test typeof(dbsplineu(4,3,5,b)[3]) == Array{Float64,2}
        @test typeof(dbsplineu(4,3,5,b)[4]) == Array{Array{Int64,1},1}
        @test isapprox(dbsplineu(4,3,5,b)[1],[1.5 3.0 4.5 6.0 7.5; 5.0 7.875 6.5 4.25 4.5; 0.0 0.0 0.0 0.0 0.0],atol=1e-5)
        @test isapprox(dbsplineu(4,3,5,b)[2],[3.0 3.0 3.0 3.0 3.0; 10.0 1.5 -7.0 -2.0 3.0; 0.0 0.0 0.0 0.0 0.0],atol=1e-5)
        @test isapprox(dbsplineu(4,3,5,b)[3],[0.0 0.0 0.0 0.0 0.0; -17.0 -17.0 10.0 10.0 10.0; 0.0 0.0 0.0 0.0 0.0],atol=1e-5)
    end

    @testset "nmatrix" begin
        @test typeof(nmatrix(3)[1]) == Array{Float64,2}
        @test typeof(nmatrix(3)[2]) == Float64
        @test nmatrix(3)[1] == [1.0 -2.0 1.0; -2.0 2.0 0.0; 1.0 1.0 0.0]
        @test nmatrix(3)[2] == 0.5
        @test nmatrix(4)[1] == [-1.0 3.0 -3.0 1.0; 3.0 -6.0 3.0 0.0; -3.0 0.0 3.0 0.0; 1.0 4.0 1.0 0.0]
        @test isapprox(nmatrix(4)[2],0.1666666666,atol=1e-10)
    end

    @testset "matpbspl" begin
        b = [0. 3 6 9; 0 10 3 6; 0 0 0 0]
        @test typeof(matpbspl(4,3,5,b)[1]) == Array{Float64,2}
        @test typeof(matpbspl(4,3,5,b)[2]) == Array{Array{Int64,1},1}
        @test isapprox(matpbspl(4,3,4,b)[1],[1.5 2.5 3.5 4.5 5.5 6.5 7.5 7.83333 6.83333 4.5 2.16667 1.16667 1.5 2.5 3.5; 5.0 7.38889 7.88889 6.5 4.72222 4.05556 4.5 5.0 4.5 3.0 1.88889 2.55556 5.0 7.38889 7.88889; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0],atol=1e-5)
    end
end