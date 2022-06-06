@time @testset "molecule.jl" begin
    @test molecule("methane.xyz") == Molecule(["C", "H", "H", "H", "H"], [[1.0, 0.3, 0.5], [1.0, 0.3, 0.5], [1.0, 0.3, 0.5], [1.0, 0.3, 0.5], [1.0, 0.3, 0.5]])
end