function θ(l, ℓᵢ, ℓⱼ, Pᵢ⃗, Pⱼ⃗, r, γ)
    θ  = cₖ(l, ℓᵢ, ℓⱼ, Pᵢ, Pⱼ)
    θ *= factorial(l) * γ^(r - l)
    θ /= factorial(r) * factorial(l - (2 * r))

    return θ
end

function gᵢ(l, lₚ, r, rₚ, i, ℓᵢ, ℓⱼ, A⃗, B⃗, P⃗, γ₁, ℓₖ, ℓᵥ, C⃗, D⃗, Q⃗, γ₂)
    δ = 1/(4 * γ₁) + 1/(4 * γ₂)

    gᵢ  = (-1)^l
    gᵢ *= θ(l, ℓᵢ, ℓⱼ, P⃗ - A⃗, P⃗ - B⃗, r, γ₁) * θ(lₚ, ℓₖ, ℓᵥ, Q⃗ - C⃗, Q⃗ - D⃗, rₚ, γ₂)
    gᵢ *= (-1)^i * (2 * δ)^(2 * (r + rₚ))
    gᵢ *= factorial(l + lₚ - (2 * r) - (2 * rₚ)) * δ^i
    gᵢ *= (P⃗ - Q⃗)^(l + lₚ - 2 * (r + rₚ + i))
    gᵢ /= (4 * δ)^(l + lₚ) * factorial(i)
    gᵢ /= factorial(l + lₚ - 2 * (r + rₚ + i))

    return gᵢ
end

function Gxyz(ℓᵢ, mᵢ, nᵢ, ℓⱼ, mⱼ, nⱼ, ℓₖ, mₖ, nₖ, ℓᵥ, mᵥ, nᵥ, i, j, k, v, Rᵢ, Rⱼ, Rₖ, Rᵥ)
    γ₁ = i + j
    γ₂ = k + v

    δ = 1/(4 * γ₁) + 1/(4 * γ₂)

    Rₚ = gaussianproduct(i, Rᵢ, j, Rⱼ, γ₁)
    Rq = gaussianproduct(k, Rₖ, v, Rᵥ, γ₂)

    IJ = distance(Rᵢ, Rⱼ)
    KV = distance(Rₖ, Rᵥ)
    PQ = distance(Rₚ, Rq)

    Gxyz = 0

    for l in 0:(ℓᵢ + ℓⱼ)
        for r in 0:trunc(Int64, (l / 2))
            for lₚ in 0:(ℓₖ + ℓᵥ)
                for rₚ in 0:trunc(Int64, (lₚ / 2))
                    for i in 0:trunc(Int64, (l + lₚ - (2 * rₚ)) / 2)
                        gx = gᵢ(l, lₚ, r, rₚ, i, ℓᵢ, ℓⱼ, Rᵢ[1], Rⱼ[1], Rₚ[1], γ₁, ℓₖ, ℓᵥ, Rₖ[1], Rᵥ[1], Rq[1], γ₂)
                    end
                end
            end
        end
    end
end