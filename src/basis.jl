abstract type Basis end
struct GaussianBasis <: Basis
    Z::Int64
    R::Int64
    l::Int64
    m::Int64
    n::Int64
    α::Vector{Float64}
end
    