function Kxyz(Rᵢ, Rⱼ, αᵢ, αⱼ, ℓᵢ, ℓⱼ, mᵢ, mⱼ, nᵢ, nⱼ)
    K = αⱼ * (2 * (ℓⱼ + mⱼ + nⱼ) + 3) * Sxyz(Rᵢ, Rⱼ, αᵢ, αⱼ, ℓᵢ, ℓⱼ, mᵢ, mⱼ, nᵢ, nⱼ)
    K -= (2 * (αⱼ^2)) * Sxyz(Rᵢ, Rⱼ, αᵢ, αⱼ, ℓᵢ, ℓⱼ + 2, mᵢ, mⱼ, nᵢ, nⱼ)
    K -= (2 * (αⱼ^2)) * Sxyz(Rᵢ, Rⱼ, αᵢ, αⱼ, ℓᵢ, ℓⱼ, mᵢ, mⱼ + 2, nᵢ, nⱼ)
    K -= (2 * (αⱼ^2)) * Sxyz(Rᵢ, Rⱼ, αᵢ, αⱼ, ℓᵢ, ℓⱼ, mᵢ, mⱼ, nᵢ, nⱼ + 2)

    K -= (1 / 2) * (ℓⱼ * (ℓⱼ - 1)) * Sxyz(Rᵢ, Rⱼ, αᵢ, αⱼ, ℓᵢ, ℓⱼ - 2, mᵢ, mⱼ, nᵢ, nⱼ)
    K -= (1 / 2) * (mⱼ * (mⱼ - 1)) * Sxyz(Rᵢ, Rⱼ, αᵢ, αⱼ, ℓᵢ, ℓⱼ, mᵢ, mⱼ - 2, nᵢ, nⱼ)
    K -= (1 / 2) * (nⱼ * (nⱼ - 1)) * Sxyz(Rᵢ, Rⱼ, αᵢ, αⱼ, ℓᵢ, ℓⱼ, mᵢ, mⱼ, nᵢ, nⱼ - 2)

    return K
end

function kinetic(basis, molecule::Molecule)
    K = length(basis)
    T = zeros(K, K)

    for (i, basisᵢ) in enumerate(basis)
        for (j, basisⱼ) in enumerate(basis)
            for (αᵢ, dᵢ) in zip(basisᵢ.α, basisᵢ.d)
                for (αⱼ, dⱼ) in zip(basisⱼ.α, basisⱼ.d)

                    Rᵢ = basisᵢ.R
                    Rⱼ = basisⱼ.R

                    ℓᵢ, mᵢ, nᵢ = basisᵢ.ℓ, basisᵢ.m, basisᵢ.n
                    ℓⱼ, mⱼ, nⱼ = basisⱼ.ℓ, basisⱼ.m, basisⱼ.n

                    T[i, j] += (
                        exp(-αᵢ * αⱼ * distance(Rᵢ, Rⱼ) / (αᵢ + αⱼ)) *
                        normalization(αᵢ, ℓᵢ, mᵢ, nᵢ) *
                        normalization(αⱼ, ℓⱼ, mⱼ, nⱼ) *
                        dᵢ *
                        dⱼ *
                        Kxyz(Rᵢ, Rⱼ, αᵢ, αⱼ, ℓᵢ, ℓⱼ, mᵢ, mⱼ, nᵢ, nⱼ)
                    )
                end
            end
        end
    end

    return T
end
