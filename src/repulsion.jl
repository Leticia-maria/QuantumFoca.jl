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

function Gxyz(ℓᵢ, mᵢ, nᵢ, ℓⱼ, mⱼ, nⱼ, ℓₖ, mₖ, nₖ, ℓᵥ, mᵥ, nᵥ, a, b, c, d, Rᵢ, Rⱼ, Rₖ, Rᵥ)
    γ₁ = a + b
    γ₂ = c + d

    δ = 1/(4 * γ₁) + 1/(4 * γ₂)

    Rₚ = gaussianproduct(a, Rᵢ, b, Rⱼ, γ₁)
    Rq = gaussianproduct(c, Rₖ, d, Rᵥ, γ₂)

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

                        for m in 0:(mᵢ + mⱼ)
                            for s in 0:trunc(Int64, (m / 2))
                                for mₚ in 0:(mₖ + mᵥ)
                                    for sₚ in 0:trunc(Int64, (mₚ / 2))
                                        for j in 0:trunc(Int64, (m + mₚ - (2 * sₚ)) / 2)
                                            gy = gᵢ(m, mₚ, s, sₚ, j, mᵢ, mⱼ, Rᵢ[2], Rⱼ[2], Rₚ[2], γ₁, mₖ, mᵥ, Rₖ[2], Rᵥ[2], Rq[2], γ₂)

                                            for n in 0:(nᵢ + nⱼ)
                                                for t in 0:trunc(Int64, (n / 2))
                                                    for nₚ in 0:(nₖ + nᵥ)
                                                        for tₚ in 0:trunc(Int64, (nₚ / 2))
                                                            for k in 0:trunc(Int64, (n + nₚ - (2 * tₚ)) / 2)
                                                                gz = gᵢ(n, nₚ, t, tₚ, k, nᵢ, nⱼ, Rᵢ[1], Rⱼ[1], Rₚ[1], γ₁, nₖ, nᵥ, Rₖ[1], Rᵥ[1], Rq[1], γ₂)

                                                                ν = l + lₚ + m + mₚ + n + nₚ - 2 * (r + rₚ + s + sₚ + t + tₚ) - (i + j + k)
                                                                F = boys(ν, (PQ / (4 * δ)))

                                                                Gxyz += gx * gy * gz * F
                                                            end
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    Gxyz *= (2 * π^2) / (γᵢ * γ₂)
    Gxyz *= sqrt(π / (γ₁ + γ₂))
    Gxyz *= exp(- (a * b * IJ) / γ₁)
    Gxyz *= exp(- (c * d * KV) / γ₂)

    Nᵢ = normalization(a, ℓᵢ, mᵢ, nᵢ)
    Nⱼ = normalization(b, ℓⱼ, mⱼ, nⱼ)
    Nₖ = normalization(c, ℓₖ, mₖ, nₖ)
    Nᵥ = normalization(d, ℓᵥ, mᵥ, nᵥ)

    Gxyz *= Nᵢ * Nⱼ * Nₖ * Nᵥ

    return Gxyz
end

function repulsion(basis. molecule)
    K = length(basis)
    V = zeros(K, K, K, K)

    Ntei = 0

    for (i, basisᵢ) in enumerate(basis)
        for (j, basisⱼ) in enumerate(basis)
            for (k, basisₖ) in enumerate(basis)
                for (v, basisᵥ) in enumerate(basis)
                    Ntei += 1

                    for (αᵢ, dᵢ) in zip(basisᵢ.α, basisᵢ.d)
                        for (αⱼ, dⱼ) in zip(basisⱼ.α, basisⱼ.d)
                            for (αₖ, dₖ) in zip(basisₖ.α, basisₖ.d)
                                for (αᵥ, dᵥ) in zip(basisᵥ.α, basisᵥ.d)

                                    Rᵢ = basisᵢ.R
                                    Rⱼ = basisⱼ.R
                                    Rₖ = basisₖ.R
                                    Rᵥ = basisᵥ.R

                                    ℓᵢ, mᵢ, nᵢ = basisᵢ.ℓ, basisᵢ.m, basisᵢ.n
                                    ℓⱼ, mⱼ, nⱼ = basisⱼ.ℓ, basisⱼ.m, basisⱼ.n
                                    ℓₖ, mₖ, nₖ = basisₖ.ℓ, basisₖ.m, basisₖ.n
                                    ℓᵥ, mᵥ, nᵥ = basisᵥ.ℓ, basisᵥ.m, basisᵥ.n

                                    tei  = dᵢ * dⱼ * dₖ * dᵥ 
                                    tei *= Gxyz(ℓᵢ, mᵢ, nᵢ, ℓⱼ, mⱼ, nⱼ, ℓₖ, mₖ, nₖ, ℓᵥ, mᵥ, nᵥ, i, j, k, v, Rᵢ, Rⱼ, Rₖ, Rᵥ)

                                    G[i, j, k, v] += tei
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    return G
end

            
                    