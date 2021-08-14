import myLibs: ComputeTasks 

using myLibs.ComputeTasks: CompTask 

module M 


Read(P;kw...) = "read"

Compute(P;kw...) = "compute" 

FoundFiles(P;kwa...) = rand(Bool)

NrParamSets = 1
	
apnt =(
										length = [10,20],
									 	width = [7],
										Barrier_height = [0,0.5],
										SCpx_magnitude = [0.6],
										)  

allparams() = Dict(k=>apnt[k] for k in keys(apnt))

end

module M2


Read(P1,P2;kwargs...) = "read2"

Compute(P1,P2; kwargs...) = "compute2" 

FoundFiles(P1, P2;kwargs...) = rand(Bool)#false#true 

NrParamSets = 2

allparams() = Dict(:x=>[1,2])
allparams(q) = Dict(:y=>[3,"4"])

end 


P = Dict(:length=>10,:width=>7, :Barrier_height=>1.0,:SCpx_magnitude=>0.4,:SCDW_position=>0.3,:delta=>0.002,:AtomToLayer=>"forced")


for (m,ps) in ([M,(P,)],[M2,(P,P,)])

	println("----------")
#	@show m ps 

	rmv(l::Int, p::AbstractDict) = Dict(k=>p[k] for k in filter(!isequal(:width),collect(keys(p)))) 

	add(l::Int, p::AbstractDict) = merge(p,Dict(:Q=>0.3))


	local task = CompTask(m; 
																rmv_internal_key = rmv,
																add_internal_param = add,
																)

#	local	get_plotparams, get_paramcombs, files_exist, get_data = task



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


cond(l::Int, p::AbstractDict) = get(p,:length,0)==10


foreach(println, task.get_paramcombs(;cond=cond))
println() 

qq = task.get_paramcombs()[1]

@show qq 


#@show methods(task.files_exist)

#@show typeof(qq) 


@show task.files_exist(qq...),task.files_exist(ps...) 

println() 

@show task.get_data(qq...)



ComputeTasks.missing_data(task; show_missing=true)

ComputeTasks.existing_data(task, show_existing=true)


ComputeTasks.get_data_one(task, mute=false)



#ComputeTasks.get_plot_one(task)


ComputeTasks.get_data_all(task, check_data=true, mute=false)

nothing 

end 















































