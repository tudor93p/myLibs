using Revise 

# test2
using Test, BenchmarkTools 


for f in [

#"mydistrib",
#"algebra" ,
#"signal",
#"identify_sectors", 
#"caroli_cond_methods",
#"taylor",
#
#"lattices_vecs",
#
#"lattices",
#
#"geometry",
#
#"lattices_attachleads",
#
#"sancho-rubio",
#
#
#"lattices_relabsdim",
#
#
#"hopping_terms",
#
#"rgf_LM08" , 
#"rgf_sq_honey_leads" ,
#"rgf_LM10" , 
#"rgf_LM11" , 
#
#"rgf_LM12",
#
#"rgf_B" , 
#
#"gf_decim_dataH",
#"layeredsystem",
#"gf_err",
#"gf",
#
#"site_transm_fromData",
#"site_transm_fromData_mem",
#
#"transm_ZHXS", 
#
#
#"transm_ZHXS_Fig2", 
#"transm_ZHXS_Fig3", 
#"transm_ZHXS_Fig5", 
#"transm_ZHXS_Fig7", 
#"transm_ZHXS_local", 
#
#
#
#"tbmodel",
#
#"TBmodel_faster-TBmatrices",
#
#"utils",
#
#"recursive_merge",
#"utils_splitws",
#"oper" ,
#
#"mirror_oper",
#
#"bandstr" ,
#
#"simple_wlos",
#
#"tasks2",
#"tasks",
#
#"tasks_multi",
#
#"param" #|> pr_in 

#"param2",
#%"rw", 
"fpos",
]

	println("\n********* $f ********* \n")

	include("$f.jl")

end









#"choose-obs" #|> pr_in



















nothing
