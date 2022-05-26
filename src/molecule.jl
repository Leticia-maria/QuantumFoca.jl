abstract type ChemicalEntity end
abstract type ChemicalFile end

struct XYZfile<:ChemicalFile
    coordsX::Matrix{Int64}
    coordsY::Matrix{Int64}
    coordsZ::Matrix{Int64}
end

struct Molecule<:ChemicalEntity
    numberOfAtoms::Int
    symbols::Vector{String}
    coords::Matrix{Int64}
end