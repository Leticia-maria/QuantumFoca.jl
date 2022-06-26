module QuantumFoca
    using LinearAlgebra
    using SpecialFunctions

    include("atoms.jl")
    include("molecule.jl")
    include("sto3g.jl")
    include("basis.jl")
    include("auxiliary.jl")
    include("overlap.jl")
    include("kinetic.jl")
    include("boys.jl")
    include("attraction.jl")
    include("repulsion.jl")
    
    export Molecule
    export molecule

    export GaussianBasis
    export buildbasis

    export doublefactorial
    export gaussianproduct
    export normalization

    export overlap
    export kinetic
    export boys
    export attraction
    export repulsion
end
