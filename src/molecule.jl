abstract type ChemicalEntity end

abstract type ChemicalFile end
struct Molecule <: ChemicalEntity
    atoms::Vector{String}
    coords::Matrix{Float64}
end

struct XYZ{M<:Molecule} <: ChemicalFile end

function Base.read(xyzfile::IO, ::Type{XYZ{M}}) where {M<:Molecule}
    elements = String[]
    coordinates = Matrix{Float64}(undef, 3, 0)

    for atom in Iterators.drop(eachline(xyzfile), 2)
        fields = split(line)

        element = fields[1]
        push!(elements, element)

        coordinate = parse.(Float64, fields[2:4])
        push!(coordinates, coordinate)
        coordinates = hcat(coordinates, coordinate)
    end

    return M(elements, coordinates)
end

Base.read(input::IO, ::Type{XYZ{Molecule}}) = read(input, XYZ{Molecule})