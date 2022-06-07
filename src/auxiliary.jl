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