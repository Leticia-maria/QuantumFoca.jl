function θ(l, lA, lB, PA, PB, r, g)
    θ  = cₖ(l, lA, lB, PA, PB)
    θ *= factorial(l) * g^(r - l)
    θ /= factorial(r) * factorial(l - 2 * r)

    return θ
end

function gi(l, lp, r, rp, i, lA, lB, Ai, Bi, Pi, gP, lC, lD, Ci, Di, Qi, gQ)
    δ = 1 / (4 * gP) + 1 / (4 * gQ)

    gi  = (-1.0)^l
    gi *= θ(l, lA, lB, Pi - Ai, Pi - Bi, r, gP) * θ(lp, lC, lD, Qi - Ci, Qi - Di, rp, gQ)
    gi *= (-1.0)^i * (2 * δ)^(2 * (r + rp))
    gi *= factorial(l + lp - 2 * r - 2 * rp) * δ^i
    gi *= (Pi - Qi)^(l + lp - 2 * (r + rp + i))
    gi /= (4 * δ)^(l + lp) * factorial(i)
    gi /= factorial(l + lp - 2 * (r + rp + i))

    return gi
end

function Gxyz(lA, mA, nA, lB, mB, nB, lC, mC, nC, lD, mD, nD, a, b, c, d, RA, RB, RC, RD)

end