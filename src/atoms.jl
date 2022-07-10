Z = Dict(
    "H" => 1,
    "He" => 2,
    "Li" => 3,
    "Be" => 4,
    "B" => 5,
    "C" => 6,
    "N" => 7,
    "O" => 8,
    "F" => 9,
    "Ne" => 10,
)

function orbitalconfig(Z::Int64)
    if Z < 3
        return ["1s"]
    else
        return ["1s", "2s", "2p"]
    end
end
