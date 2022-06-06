@time @testset "molecule.jl" begin
    methane = molecule("methane.xyz")
    @test methane.atoms == ["C", "H", "H", "H", "H"]
end