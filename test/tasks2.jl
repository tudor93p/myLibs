import myLibs: ComputeTasks, Parameters 

using myLibs.ComputeTasks: CompTask 
using myLibs.Parameters: Calculation

using OrderedCollections: OrderedDict 
using Combinatorics: powerset 



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#
input_dict = Dict(

		:allparams => Dict(:X=>[1,8],:Y=>[3,4],:T=>[10,20,30]),

#		:digits => OrderedDict(:X=>(1,0), :Y=>(1,1), :Z=>(2,2)),
		:digits => (X=(1,0), Y=(1,1), Z=(2,2)),

		)




module M 

import ..input_dict,..Parameters

Read(P;kw...) = "read"

Compute(P;kw...) = [sleep(rand()*2);"compute"]

FoundFiles(P;kwa...) = rand(Bool)

NrParamSets = 1


	#export usedkeys
	usedkeys(args...) = [:X]
	a=3 
	
#	NrParamSets = 2 
	
	digits = (X=(3,3),) 
	
	path = "abc" 
	
	params_digits = Parameters.typical_params_digits(usedkeys, input_dict[:digits])
	
	allparams = Parameters.typical_allparams(usedkeys, input_dict[:allparams])

	
end  


PF = Parameters.ParamFlow(M.path, M, M.allparams)




C = Calculation(PF, M)
	
task = CompTask(C)
 

ps = [Dict(:length=>10,:width=>7, :Barrier_height=>1.0,:SCpx_magnitude=>0.4,:SCDW_position=>0.3,:delta=>0.002,:AtomToLayer=>"forced")]


for (k,v) in task.get_plotparams(ps...)


	println(k, " ", v)

end 

@show task.get_plotparams()




println() 


for pc in task.get_paramcombs()

	foreach(println, pc)

	println()

end

println()  



qq = task.get_paramcombs()[1]

@show qq 




@show task.files_exist(qq...),task.files_exist(ps...) 

println() 

@show task.get_data(qq...)



ComputeTasks.missing_data(task; show_missing=true)

ComputeTasks.existing_data(task, show_existing=true)


ComputeTasks.get_data_one(task, mute=false)



#ComputeTasks.get_plot_one(task)


ComputeTasks.get_data_all(task, check_data=false,  mute=true)

nothing 


