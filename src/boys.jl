function boys(ν, x)
    if x < 1e-7
        return (2 * ν + 1)^(-1.0) - x * (2 * ν + 3)^(-1.0)
    else
        return (1 / 2) * x^(-(ν + 0.5)) * gamma(ν + 0.5) * gamma_inc(ν + 0.5, x)[1]
    end
end
