abstract type Basis end

"""
```GaussianBasis``` is a *subtype* of ```Basis``` that stores the coefficients, exponents and angular momenta of the atomic orbital.
A basis set in theoretical and computational chemistry is a set of functions (called basis functions) that is used to represent the 
electronic wave function in the Hartree–Fock method or density-functional theory in order to turn the partial differential equations 
of the model into algebraic equations suitable for efficient implementation on a computer.
The use of basis sets is equivalent to the use of an approximate resolution of the identity: 
the orbitals ``|\psi _{i}\rangle`` are expanded within the basis set as a linear combination 
of the basis functions ``|\psi _{i}\rangle \approx \sum _{\mu }c_{\mu i}|\mu \rangle``, where the expansion coefficients 
``c_{\mu i}`` are given by ``c_{\mu i}=\sum _{\nu }\langle \mu |\nu \rangle ^{-1}\langle \nu |\psi _{i}\rangle``.
The basis set can either be composed of atomic orbitals (yielding the linear combination of atomic orbitals approach), which is the 
usual choice within the quantum chemistry community; plane waves which are typically used within the solid state community, or real-space 
approaches. Several types of atomic orbitals can be used: Gaussian-type orbitals, Slater-type orbitals, or numerical atomic orbitals.
Out of the three, Gaussian-type orbitals are by far the most often used, as they allow efficient implementations of Post-Hartree–Fock methods. 
"""
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