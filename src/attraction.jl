function vᵢ(l, r, i, ℓᵢ, ℓⱼ, Aᵢ, Bᵢ, Cᵢ, Pᵢ, γ)
    ϵ = 1 / (4 * γ)

    vᵢ  = (-1)^l
    vᵢ *= cₖ(l, ℓᵢ, ℓⱼ, Pᵢ - Aᵢ, Pᵢ - Bᵢ)
    vᵢ *= (-1)^i * factorial(l)
    vᵢ *= (Pᵢ - Cᵢ)^(l - (2 * r) - (2 * i)) * ϵ^(r + i)
    vᵢ /= factorial(r)
    vᵢ /= factorial(i)
    vᵢ /= factorial(l - (2 * r) - (2 * i))

    return vᵢ
end

function Vxyz(ℓᵢ, mᵢ, nᵢ, ℓⱼ, mⱼ, nⱼ, αᵢ, αⱼ, Rᵢ, Rⱼ, Rₖ, Z)
    γ = αᵢ + αⱼ

    Rₚ = gaussianproduct(αᵢ, Rᵢ, αⱼ, Rⱼ, γ)

    IJ = distance(Rᵢ, Rⱼ)
    PK = distance(Rₚ, Rₖ)

    Vxyz = 0
    
    for l in 0:(ℓᵢ + ℓⱼ)
        for r in 0:trunc(Int64, (l / 2))
            for i in 0:trunc(Int64, ((l - (2 * r)) / 2))
                Vx = vᵢ(l, r, i, ℓᵢ, ℓⱼ, Rᵢ[1], Rⱼ[1], Rₖ[1], Rₚ[1], γ)

                for m in 0:(mᵢ + mⱼ)
                    for s in 0:trunc(Int64, (m / 2))
                        for j in 0:trunc(Int64, ((m - (2 * s)) / 2))
                            Vy = vᵢ(m, s, j, mᵢ, mⱼ, Rᵢ[2], Rⱼ[2], Rₖ[2], Rₚ[2], γ)

                            for n in 0:(nᵢ + nⱼ)
                                for t in 0:trunc(Int64, (n / 2))
                                    for k in 0:trunc(Int64, ((n - (2 * t)) / 2))
                                        Vz = vᵢ(n, t, k, nᵢ, nⱼ, Rᵢ[3], Rⱼ[3], Rₖ[3], Rₚ[3], γ)

                                        ν = l + m + n - 2 * (r + s + t) - (i + j + k)
                                        F = boys(ν, (γ * abs(PK)))

                                        Vxyz += Vx * Vy * Vz * F
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    Nᵢ = normalization(αᵢ, ℓᵢ, mᵢ, nᵢ)
    Nⱼ = normalization(αⱼ, ℓⱼ, mⱼ, nⱼ)
    
    Vxyz *= (2 * π) / γ
    Vxyz *= exp(-αᵢ * αⱼ * abs(IJ) / γ)
    Vxyz *= Nᵢ * Nⱼ
    Vxyz *= -Z

    return Vxyz
end

function attraction(basis, molecule::Molecule)
    n = length(molecule.atoms)
    K = length(basis)
    V = zeros(K, K, n)

    for (i, basisᵢ) in enumerate(basis)
        for (j, basisⱼ) in enumerate(basis)
            for (k, Z) in enumerate(molecule.numbers)
                for (αᵢ, dᵢ) in zip(basisᵢ.α, basisᵢ.d)
                    for (αⱼ, dⱼ) in zip(basisⱼ.α, basisⱼ.d)

                        Rᵢ = basisᵢ.R
                        Rⱼ = basisⱼ.R
                        Rₖ = molecule.coords[k, :]
    
                        ℓᵢ, mᵢ, nᵢ = basisᵢ.ℓ, basisᵢ.m, basisᵢ.n
                        ℓⱼ, mⱼ, nⱼ = basisⱼ.ℓ, basisⱼ.m, basisⱼ.n
                        
                        V[i, j, k] += dᵢ * dⱼ *
                                      Vxyz(ℓᵢ, mᵢ, nᵢ, ℓⱼ, mⱼ, nⱼ, αᵢ, αⱼ, Rᵢ, Rⱼ, Rₖ, Z)
                    end
                end
            end
        end
    end

    return V
end