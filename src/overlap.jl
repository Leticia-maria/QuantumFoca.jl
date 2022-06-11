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

    for i in 1:K, j in 1:K
        ao₁ = basis[i]
        ao₂ = basis[j]
        S[i, j] += (exp.(-ao₁.α .* ao₂.α .* 
                    distance(ao₁.R, ao₂.R) ./ 
                    (ao₁.α .+ ao₂.α)) .*
                    normalization.(ao₁.α, ao₁.ℓ, ao₁.m, ao₁.n) * 
                    normalization.(ao₂.b, ao₂.ℓ, ao₂.m, ao₂.n) *
                    ao₁.d .* ao₂.d .* 
                    Sxyz(ao₁.R, ao₂.R, ao₁.α, ao₂.α, ao₁.ℓ, ao₂.ℓ, ao₁.m, ao₂.m, ao₁.n, ao₂.n))
    end
    return S
end