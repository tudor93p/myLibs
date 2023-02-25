using Revise 

# test2
using Test, BenchmarkTools 


for f in [

#"mydistrib",
#"algebra" ,
 
"signal",
#
#
#"identify_sectors", 
#
#
#"caroli_cond_methods"


#"taylor",

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
#"gf_err",
#"gf",


#"transm_ZHXS", 


#"transm_ZHXS_Fig2", 
#"transm_ZHXS_Fig3", 
#"transm_ZHXS_Fig5", 
#"transm_ZHXS_Fig7", 
#"transm_ZHXS_local", 



#"tbmodel",

#"utils",

#"utils_splitws",
#"oper" ,

#"mirror_oper",

#"bandstr" |> pr_in



#"tasks_multi",

##"param" #|> pr_in 

#"param2",

]

	println("\n********* $f ********* \n")

	include("$f.jl")

end





#"tasks" #|> pr_in




#"choose-obs" #|> pr_in



















nothing
