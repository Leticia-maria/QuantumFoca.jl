abstract type Basis end
struct GaussianBasis<:Basis
    atomNo::Int
    orbital::String
    coords::Array{Float, 3}
    l::Int
    m::Int
    n::Int
    a::Array{Float, 3}
    d::Array{Float, 3}
end 

function build_sto3Gbasis(molecule::Molecule)
    sto3Gbasis = []
    K = 0

    
end