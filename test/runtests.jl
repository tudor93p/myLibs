using Revise 

# test2
using Test  


for f in (

#"algebra" |> pr_in

#"lattices_vecs" |> pr_in

#"lattices",

#"geometry",

#"lattices_attachleads",

#"sancho-rubio" |> pr_in 


#"lattices_relabsdim" |> pr_in


#"hopping_terms",

#"rgf_LM08" |> pr_in 
#"rgf_sq_honey_leads" |> pr_in
#"rgf_LM10" |> pr_in 
#"rgf_LM11" |> pr_in 

#"rgf_LM12",

#"rgf_B" |> pr_in 

#"layeredsystem",
#"gf",



#"transm_ZHXS", 

#"transm_ZHXS_Fig2", 
#"transm_ZHXS_Fig3", 
#"transm_ZHXS_Fig5", 
#"transm_ZHXS_Fig7", 
#"transm_ZHXS_local", 



#"tbmodel" |> pr_in 

#"utils",

#"utils_splitws",
"oper" ,

#"bandstr" |> pr_in




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
