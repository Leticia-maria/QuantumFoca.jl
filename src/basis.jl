abstract type Basis end
struct GaussianBasis <: Basis
    Z::Int64
    R
    l
    m
    n
    α
    