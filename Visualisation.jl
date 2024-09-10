valid_input = false

while !valid_input
    println("Do you want to obtain animation with [s]tatic streamlines, [d]ynamic streamlines or [b]oth?")
    answer = readline()

    # case-insensitive comparison
    answer = lowercase(answer)

    if answer == "s"
        println("You have selected static streamlines.")
        valid_input = true
        include("visualisation_Schwarzschild.jl")
    elseif answer == "d"
        println("You have selected dynamic streamlines.")
        valid_input = true
        include("visualisation_Schwarzschild_vanishing.jl")
    elseif answer == "b"
        println("You have selected both static and dynamic streamlines.")
        valid_input = true
        include("visualisation_Schwarzschild.jl")
        include("visualisation_Schwarzschild_vanishing.jl")
    else
        println("Invalid input. Please enter 's', 'd', or 'b'.")
    end
end


