import myLibs: ComputeTasks 


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

	@show m ps 

	rmv(l,p::AbstractDict) = Dict(k=>p[k] for k in filter(!isequal(:width),collect(keys(p))))
	add(l,p::AbstractDict) = merge(p,Dict(:Q=>0.3))


	local task = ComputeTasks.init_task(M; 
																rmv_internal_key = rmv,
																add_internal_param = add,
																)

	local	get_plotparams, get_paramcombs, files_exist, get_data = task



for (k,v) in get_plotparams(ps...)

	println(k, " ", v)

end 

println() 



foreach(println, get_paramcombs())

println()  


cond(p::AbstractDict,l) = get(p,:length,0)==10


foreach(println, get_paramcombs(;cond=cond))
println() 

@show 	files_exist(get_paramcombs()[1]...)
@show get_data(get_paramcombs()[1]...)


end 















































nothing
