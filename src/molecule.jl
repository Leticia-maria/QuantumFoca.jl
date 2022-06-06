abstract type ChemicalEntity end
struct Molecule <: ChemicalEntity
    atoms::Vector{String}
    coords::Matrix{Float64}
end

function molecule(xyzfile::String) :: Molecule
    elements = String[]
    coordinates = []

    for line in Iterators.drop(eachline(xyzfile), 2)
        fields = split(line)

        element = fields[1]
        push!(elements, element)

        coordinate = parse.(Float64, fields[2:4])
        push!(coordinates, coordinate)
    end
    
    coordinates = mapreduce(permutedims, vcat, coordinates)
    return Molecule(elements, coordinates)
end