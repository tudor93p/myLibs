using Revise 

# test2
using Test 

function pr_in(f::AbstractString) 

	println("\n********* $f ********* \n")

	include("$f.jl")

end
#"algebra" |> pr_in

#"lattices_vecs" |> pr_in

#"lattices" |> pr_in 

#"geometry" |> pr_in 

#"lattices_attachleads" |> pr_in

#"sancho-rubio" |> pr_in 


#"lattices_relabsdim" |> pr_in

"rgf_LM08" |> pr_in 
#"rgf_sq_honey_leads" |> pr_in
#"rgf_LM10" |> pr_in 
#"rgf_LM11" |> pr_in 

#"rgf_LM12" |> pr_in 

#"rgf_B" |> pr_in 


#"rgf" |> pr_in
#"gf" |> pr_in


#"tbmodel" |> pr_in 

#"utils" |> pr_in 

#"utils_splitws" |> pr_in 
#"oper" |> pr_in

#"bandstr" |> pr_in



#"layeredsystem" |> pr_in

##"param" #|> pr_in 

#"param2" |>pr_in 






#"tasks" #|> pr_in


#"tasks_multi" |> pr_in


#"choose-obs" #|> pr_in



















nothing
