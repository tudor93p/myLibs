using Revise 

using Test 


function pr_in(f::AbstractString) 

	println("\n********* $f ********* \n")

	include("$f.jl")

end
#include("algebra.jl")


#"lattices" |> pr_in 

#"geometry" |> pr_in 

#"lattices_attachleads" |> pr_in

#
#include("tbmodel.jl")
#
#include("utils.jl")
#

#include("bandstr.jl") 
#
#

#include("layeredsystem.jl")

#"param" |> pr_in 
#
"param2" |>pr_in
#
#
#
#
#"tasks" |> pr_in
#
#
#
#
#
#
#















nothing
