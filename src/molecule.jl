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

"""
Build 2D nested list containing raw data from .xyz input file.
"""
function constructNestedArray(data)
    n = 0
    fileNestedlist = Matrix{Int64}
    open(data) do file
        for (l, line) in enumerate(eachline(file))
            if l == 1
                try
                    n = parse(Int64, line)
                catch
                    print("The xyz file is not in the correct format. Make sure the format follows: https://en.wikipedia.org/wiki/XYZ_file_format")
                end
            end

            if l  == 2
                continue
            end

            if l != 1
                atom_coordinates = [
                    split(line)[1],
                    parse(Float64, split(line)[2]),
                    parse(Float64, split(line)[3]),
                    parse(Float64, split(line)[4])
                ]
                append!(fileNestedlist, [atom_coordinates])
            end
        end
    end
    return n, fileNestedlist
end