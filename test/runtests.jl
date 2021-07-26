using Revise 

using Test 

#include("algebra.jl")
#include("lattices.jl")
#
#include("tbmodel.jl")
#
#include("utils.jl")
#

#include("bandstr.jl") 
#
#

function pr_in(f::AbstractString) 

	println("\n********* $f ********* \n")

	include("$f.jl")

end

#include("layeredsystem.jl")

"param" |> pr_in 

"param2" |>pr_in




"tasks" |> pr_in






















nothing
