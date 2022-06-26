push!(LOAD_PATH,"../src/")

using QuantumFoca
using Documenter

makedocs(
         sitename = "QuantumFoca.jl",
         modules  = [QuantumFoca],
         pages=[
                "Home" => "index.md",
                "API" => [
                    "Input" => "input.md",
                    "Basis Sets" => "basis.md"
                ]
               ])
               
deploydocs(;
    repo="github.com/Leticia-maria/QuantumFoca.jl.git"
)