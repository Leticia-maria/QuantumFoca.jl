push!(LOAD_PATH,"../src/")

using Foca
using Documenter

makedocs(
         sitename = "Foca.jl",
         modules  = [Foca],
         pages=[
                "Home" => "index.md",
                "API" => [
                    "Input" => "input.md",
                    "Basis Sets" => "basis.md"
                ]
               ])
               
deploydocs(;
    repo="github.com/Leticia-maria/Foca.jl.git"
)