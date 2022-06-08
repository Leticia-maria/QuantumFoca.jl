"""
This function calculate the double factorial of a number.
"""
function doublefactorial(number)
    fact = one(number)

    if number%2==0
        for m in 1:number
            if m%2 == 0
                fact *= m
            end
        end
    elseif number%2==1
        for m in 1:number
            if m%2== 1
                fact *= m
            end
        end
    end
    return fact
end

"""
gaussian_norm(α, l, m, n)
Compute the normalization constant for a Gaussian primitive basis function.
"""
@inline function gaussian_norm(α, l::Int, m::Int, n::Int)
    N = sqrt((2α / π)^3) * (4α)^(l + m + n)
    κ = doublefactorial(2l - 1) * doublefactorial(2m - 1) * doublefactorial(2n - 1)

    return sqrt(N / κ)
end