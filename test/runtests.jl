using Revise 

# test2
using Test  


for f in (

#"algebra" |> pr_in

#"lattices_vecs" |> pr_in

#"lattices" |> pr_in 

#"geometry" |> pr_in 

#"lattices_attachleads" |> pr_in

#"sancho-rubio" |> pr_in 


#"lattices_relabsdim" |> pr_in


#"hopping_terms",

#"rgf_LM08" |> pr_in 
#"rgf_sq_honey_leads" |> pr_in
#"rgf_LM10" |> pr_in 
#"rgf_LM11" |> pr_in 

#"rgf_LM12",

#"rgf_B" |> pr_in 


#"transm_ZHXS", 

"transm_ZHXS_Fig2", 

#"rgf" |> pr_in
#"gf",


#"tbmodel" |> pr_in 

#"utils" |> pr_in 

#"utils_splitws" |> pr_in 
#"oper" ,

#"bandstr" |> pr_in



#"layeredsystem" |> pr_in

##"param" #|> pr_in 

#"param2" |>pr_in 

)


	println("\n********* $f ********* \n")

	include("$f.jl")

end





#"tasks" #|> pr_in


#"tasks_multi" |> pr_in


#"choose-obs" #|> pr_in



















nothing
