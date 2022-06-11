module Foca
    include("atoms.jl")
    include("molecule.jl")
    include("sto3g.jl")
    include("basis.jl")
    include("auxiliary.jl")
    include("overlap.jl")
    
    export Molecule
    export molecule

    export GaussianBasis
    export buildbasis

    export doublefactorial
    export gaussianproduct
    export normalization
    export overlap
end
