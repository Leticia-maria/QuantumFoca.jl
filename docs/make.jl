using Foca
using Documenter

makedocs(
         sitename = "Foca.jl",
         modules  = [Foca],
         pages=[
                "Home" => "src/index.md"
               ])
               
deploydocs(;
    repo="github.com/USERNAME/Foca.jl",
)