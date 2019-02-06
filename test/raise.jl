using NURBS
using Base.Test

@testset "Order Increasing" begin
    
    @testset "From order 2 to 3" begin
        x = knot(4,2)
        b = [1. 2 3 4;1 2 2 1;1 1 1 1]
        @test typeof(raise23(b,x,4)) == Tuple{Array{Float64,2},Int64}
        @test raise23(b,x,4)[1] == [1 1.5 2 2.5 3 3.5 4; 1 1.5 2 2 2 1.5 1; 1 1 1 1 1 1 1]
        @test raise23(b,x,4)[2] == 7
        x = [0, 0., 1, 1, 2, 2]
        @test typeof(raise23(b,x,4)) == Tuple{Array{Float64,2},Int64}
        @test raise23(b,x,4)[1] == [1 1.5 2 3 3.5 4; 1 1.5 2 2 1.5 1; 1 1 1 1 1 1]
        @test raise23(b,x,4)[2] == 6
    end
    
    @testset "From order 4 to 5" begin
        x = knot(4,4)
        b = [1. 2 3 4 ; 1 2 2 1 ; 1 1 1 1]
        @test typeof(raise45(b,x,4)) == Tuple{Array{Float64,2},Int64}
        @test raise45(b,x,4)[1] == [1 1.75 2.5 3.25 4; 1 1.75 2 1.75 1; 1 1 1 1 1]
        @test raise45(b,x,4)[2] == 5
        x = knot(5,4)
        b = [1. 2 3 4 5; 1 2 2 1 0; 1 1 1 1 1]
        @test typeof(raise45(b,x,5)) == Tuple{Array{Float64,2},Int64}
        @test raise45(b,x,5)[1] == [1 1.75 2.25 3 3.75 4.25 5; 1 1.75 2 1.875 1.25 0.75 0.0; 1 1 1 1 1 1 1]
        @test raise45(b,x,5)[2] == 7
        x = knot(6,4)
        b = [1. 2 3 4 5 6; 1 2 2 1 0 0; 1 1 1 1 1 1]
        @test typeof(raise45(b,x,6)) == Tuple{Array{Float64,2},Int64}
        @test isapprox(raise45(b,x,6)[1], [1 1.75 2.25 2.95833 3.5 4.04167 4.75 5.25 6; 1 1.75 2 1.91667 1.5 0.958333 0.25 0.0 0.0; 1 1 1 1 1 1 1 1 1], atol=5e-4)
        @test raise45(b,x,6)[2] == 9
        x = knot(7,4)
        b = [1. 2 3 4 5 6 7; 1 2 2 1 0 0 1; 1 1 1 1 1 1 1]
        @test typeof(raise45(b,x,7)) == Tuple{Array{Float64,2},Int64}
        @test isapprox(raise45(b,x,7)[1], [1 1.75 2.25 2.95833 3.5 4 4.5 5.04167 5.75 6.25 7; 1 1.75 2 1.91667 1.5 1 0.5 0.08333 0 0.25 1; 1 1 1 1 1 1 1 1 1 1 1]; atol=1e-5)
        @test raise45(b,x,7)[2] == 11
    end
end
