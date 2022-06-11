abstract type Basis end
struct GaussianBasis <: Basis
    R::Vector{Float64}
    α::Vector{Float64}
    d::Vector{Float64}
    ℓ::Int64
    m::Int64
    n::Int64
end

function buildbasis(molecule::Molecule)
    sto3g = GaussianBasis[]

    for i in eachindex(molecule.atoms)
        number = molecule.numbers[i]
        coord = molecule.coords[i, :]

        for orbital in orbitalconfig(number)
            if orbital == "1s"
                push!(sto3g, GaussianBasis(
                    coord,
                    sto3g_α[number][1, :],
                    sto3g_d[1, :],
                    0,
                    0,
                    0
                ))
            elseif orbital == "2s"
                push!(sto3g, GaussianBasis(
                    coord,
                    sto3g_α[number][2, :],
                    sto3g_d[2, :],
                    0,
                    0,
                    0
                ))
            elseif orbital == "2p"
                push!(sto3g, GaussianBasis(
                    coord,
                    sto3g_α[number][2, :],
                    sto3g_d[3, :],
                    1,
                    0,
                    0
                ))
                push!(sto3g, GaussianBasis(
                    coord,
                    sto3g_α[number][2, :],
                    sto3g_d[3, :],
                    0,
                    1,
                    0
                ))
                push!(sto3g, GaussianBasis(
                    coord,
                    sto3g_α[number][2, :],
                    sto3g_d[3, :],
                    0,
                    0,
                    1
                ))
            end
        end
    end

    return sto3g
end