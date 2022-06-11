using Foca
using Test

@testset "Foca.jl" begin
    include("molecule.jl")
    include("basis.jl")
    include("doublefactorial.jl")
    include("normalization.jl")
end
