using Revise 

# test2
using Test 

function pr_in(f::AbstractString) 

	println("\n********* $f ********* \n")

	include("$f.jl")

end
#include("algebra.jl")

#"lattices_vecs" |> pr_in

#"lattices" |> pr_in 

#"geometry" |> pr_in 

#"lattices_attachleads" |> pr_in

#
#"tbmodel" |> pr_in 

#"utils" #|> pr_in 

#"utils_splitws" |> pr_in 

#include("bandstr.jl") 



#"layeredsystem" |> pr_in

##"param" #|> pr_in 

#"param2" |>pr_in




#"tasks" #|> pr_in


"tasks_multi" |> pr_in


#"choose-obs" #|> pr_in



















nothing
