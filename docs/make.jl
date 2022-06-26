push!(LOAD_PATH,"../src/")

using Foca
using Documenter

makedocs(
         sitename = "Foca.jl",
         modules  = [Foca],
         pages=[
                "Home" => "index.md"
               ])
               
deploydocs(;
    repo="github.com/USERNAME/Foca.jl",
)