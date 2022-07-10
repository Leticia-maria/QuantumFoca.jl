@time @testset "auxiliary.jl" begin
    @test doublefactorial(1) == 1
    @test doublefactorial(2) == 2
    @test doublefactorial(3) == 3
    @test doublefactorial(4) == 8
    @test doublefactorial(5) == 15
    @test doublefactorial(9) == 945
end
