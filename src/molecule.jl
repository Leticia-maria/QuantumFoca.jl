""" 
```ChemicalEntity``` is an *abstract type* that englobes ```Molecule```, ```Atom```, ```AtomicOrbital``` and ```MolecularOrbital``` structures.
"""
abstract type ChemicalEntity end

"""
```Molecule``` is a *subtype* of ```ChemicalEntity``` that stores coordinates, atomic symbols and atomic numbers as objects.
"""
struct Molecule <: ChemicalEntity
    atoms::Vector{String}
    coords::Matrix{Float64}
    numbers::Vector{Int64}
end


"""
This method takes an ```.xyz``` file (with cartesian coordinates of atoms in a molecule) and returns a ```Molecule```. The ```.xyz``` file should be formatted as follows

```julia
2 

H      -1.788131055      0.000000000     -4.028513155
H      -1.331928651      0.434077746     -3.639854078
```

In the first line, the file should contain the numer of atoms that are in the molecule. In the second line, there is a comment, which can be the *name of the compound*, 
*molecular formula*, etc. To further information about ```.xyz``` files, [*click here*](https://www.reviversoft.com/file-extensions/xyz). For example, if we take the 
example file ```h2.xyz```, it is possible to give it as an input by calling ```molecule``` method.

```julia
molecule("h2.xyz")
```

The example above works if the file is in the current directory that you are working on. In other case, you can just give the path to the file of interest.

```julia
molecule(PATH)
```
"""
function molecule(xyzfile::String) :: Molecule
    elements = String[]
    coordinates = []
    Zvalues = Int[]

    for line in Iterators.drop(eachline(xyzfile), 2)
        fields = split(line)

        element = fields[1]
        push!(elements, element)
        push!(Zvalues, Z[element])

        coordinate = parse.(Float64, fields[2:4])
        push!(coordinates, coordinate)
    end
    
    coordinates = mapreduce(permutedims, vcat, coordinates)

    return Molecule(elements, coordinates, Zvalues)
end    