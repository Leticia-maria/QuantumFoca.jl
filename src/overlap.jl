function Sxyz(Rᵢ, Rⱼ, αᵢ, αⱼ, ℓᵢ, ℓⱼ, mᵢ, mⱼ, nᵢ, nⱼ)
    Rₚ = gaussianproduct(αᵢ, Rᵢ, αⱼ, Rⱼ, αᵢ + αⱼ)

    Sx = sᵢ(ℓᵢ, ℓⱼ, αᵢ + αⱼ, Rᵢ[1], Rⱼ[1], Rₚ[1])
    Sy = sᵢ(mᵢ, mⱼ, αᵢ + αⱼ, Rᵢ[2], Rⱼ[2], Rₚ[2])
    Sz = sᵢ(nᵢ, nⱼ, αᵢ + αⱼ, Rᵢ[3], Rⱼ[3], Rₚ[3])
    
    return Sx*Sy*Sz
end

function overlap(basis, molecule::Molecule)
    K = length(basis)
    S = zeros(K, K)

    for (i, basisᵢ) in enumerate(basis)
        for (j, basisⱼ) in enumerate(basis)
            for (αᵢ, dᵢ) in zip(basisᵢ.α, basisᵢ.d)
                for (αⱼ, dⱼ) in zip(basisⱼ.α, basisⱼ.d)

                    Rᵢ = basisᵢ.R
                    Rⱼ = basisⱼ.R

                    ℓᵢ, mᵢ, nᵢ = basisᵢ.ℓ, basisᵢ.m, basisᵢ.n
                    ℓⱼ, mⱼ, nⱼ = basisⱼ.ℓ, basisⱼ.m, basisⱼ.n

                    S[i, j] += (exp(-αᵢ * αⱼ * distance(Rᵢ, Rⱼ) / (αᵢ + αⱼ)) *
                                normalization(αᵢ, ℓᵢ, mᵢ, nᵢ) *
                                normalization(αⱼ, ℓⱼ, mⱼ, nⱼ) *
                                dᵢ * dⱼ *
                                Sxyz(Rᵢ, Rⱼ, αᵢ, αⱼ, ℓᵢ, ℓⱼ, mᵢ, mⱼ, nᵢ, nⱼ))
                end
            end
        end
    end
    
    return S
end