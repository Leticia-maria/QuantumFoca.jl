module Foca
    include("atoms.jl")
    include("molecule.jl")
    include("sto3g.jl")
    include("basis.jl")
    
    export Molecule
    export molecule

    export GaussianBasis
    export buildbasis
end
