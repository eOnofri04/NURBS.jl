using NURBS
using Base.Test

@testset "Base Generator" begin
    
    @testset "basis" begin
        x = knot(4,4)
        @test typeof(basis(4,4,0.5,x)) == Array{Float64,1}
        @test isapprox(basis(4,4,0.5,x),[0.125, 0.375, 0.375, 0.125],atol=1e-5)
    end
    
    @testset "dbasis" begin
        x = knot(4,4)
        @test typeof(dbasis(4,4,0.5,x)) == Tuple{Array{Float64,1}, Array{Float64,1}, Array{Float64,1}}
        @test isapprox(dbasis(4,4,0.5,x)[1],[0.125, 0.375, 0.375, 0.125],atol=1e-5)
        @test isapprox(dbasis(4,4,0.5,x)[2],[-0.75, -0.75, 0.75, 0.75],atol =1e-5)
        @test isapprox(dbasis(4,4,0.5,x)[3],[3.0, -3.0, -3.0, 3.0],atol =1e-5)
    end
    
    @testset "dbasisu" begin
        x = knotu(5,3)
        @test typeof(dbasisu(5,3,0.5,x)) == Tuple{Array{Float64,1}, Array{Float64,1}, Array{Float64,1}}
        @test isapprox(dbasisu(5,3,0.5,x)[1],[0.125, 0.0, 0.0, 0.0, 0.0],atol=1e-5)
        @test isapprox(dbasisu(5,3,0.5,x)[2],[0.5, 0.0, 0.0, 0.0, 0.0],atol=1e-5)
        @test isapprox(dbasisu(5,3,0.5,x)[3],[1.0, 0.0, 0.0, 0.0, 0.0],atol=1e-5)
    end
end