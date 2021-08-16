module ComputeTasks
#############################################################################

using Distributed

import ..Random, PyCall, ..OrderedDict 

import ..Utils, ..Parameters

using ..Parameters: UODict 



export CompTask

#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#






#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

struct CompTask

	name::String

	get_plotparams::Function 
	
	get_paramcombs::Function 
	
	files_exist::Function 
	
	get_data::Function 
	
end 


function CompTask(args::Vararg{<:Function,4})

	CompTask("", args...)

end 


get_taskname(M::Module) = string(Base.fullname(M)[end]) 
get_taskname(t) = t.name

#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#





function CompTask(M;

									found_files = M.FoundFiles,

								  read_data = M.Read,

					 			 	calc_data = M.Compute,

									valid_pcomb = nothing,
					
									kwargs...

									)::CompTask

	#found_files, read_data, calc_data take as arguments one or more 'UODict'. 
	# To be used as f(P::UODicts) or f(P...)

	rmv_internal_key, add_internal_param = Parameters.rmv_add_summarize(; kwargs...)


	get_plotparams = Parameters.f_get_plotparams(M, rmv_internal_key)

	fill_internal = Parameters.f_fill_internal(M, add_internal_param)





	function calc_or_read(p::Utils.List, force_comp::Bool, mute::Bool; 
												kwargs...)

		F,text = if !force_comp && found_files(p...; kwargs...)
			
					(read_data,"Reading files")

				else 

					(calc_data,"Calculating")

				end  


		mute && return F(p...; kwargs...)

		println()

		@info string(text, " ", Parameters.tostr(M))
		
		println()
		
	

		for d in p 
			
			Utils.isList(d) ? println.(Utils.NT.(d)) : println(Utils.NT(d))

		end 

	
		out = @time F(p...; kwargs...)
	
	  println("Finished.\n")
	
	  return out
	
	end
	



	return CompTask(

		get_taskname(M), 

		get_plotparams, 
	
	
		function get_paramcombs(;cond=nothing, repl=nothing)
	
			Parameters.get_paramcombs(M;
																
				repl = (rmv_internal_key, repl),
	
				cond = (valid_pcomb, cond),
				
				)
		end,
	

	
		function files_exist(Ps...; target=nothing)

				found_files(fill_internal(Ps)...; target=target) 
			
		end,
	
	
	
		function get_data(P...; target=nothing, 
											force_comp::Bool=false,
											mute::Bool=true, 
											shuffle::Bool=false, 
											fromPlot::Bool=false,
											get_good_P::Bool=false, 
											apply_rightaway::Function=(d,p)->d)
	
			good_P = fill_internal(
									fromPlot ? Parameters.convertParams_fromPlot(M,P...) : P
									)
	
			data = if target=="None" nothing else 
	
					apply_rightaway(calc_or_read(good_P, force_comp, mute;
																			 target=target),
													good_P)
	
			end
	
	
			return get_good_P ? (data, good_P) : data 
	
		end, 


	

 )

end

					






#===========================================================================#
#
#	Compute for all parameter combinations 
#							(optional: for which there is no data already stored)
#
#---------------------------------------------------------------------------#


function get_data_all(task; shuffle=false, seed=nothing, rev=false,
														check_data=true, mute=true) 


	AC = task.get_paramcombs()

	check_data && filter!(c->!task.files_exist(c...), AC)

	rev && reverse!(AC)


  if shuffle
	
		seed isa Int && Random.seed!(seed)

		Random.shuffle!(AC)

	end

  Utils.Distribute_Work(AC, task.get_data;
												vararg = true,
												force_comp =! check_data, 
												mute=mute)

	return

end



#===========================================================================#
#
#	Gives the number of parameters that are yet to be computed.
#			Optionally, these are printed explicitly
#
#---------------------------------------------------------------------------#

function missing_data(task; show_missing=false)

	allcombs = task.get_paramcombs()

	id = show_missing ? isequal(first(workers())) : x->false

	notdonecombs = allcombs[
													
		(nprocs()>1 ? pmap : map)(enumerate(allcombs)) do (i,c)


			out = !task.files_exist(c...)


#			id(myid()) && println(i,"/",length(allcombs)," Files ", out ? "don't " : "", "exist")

			return out

		end ]




	println(string("\nTask: ",get_taskname(task),
								 "\nTotal number of jobs left to do: ",
									length(notdonecombs),"/",length(allcombs),"\n"))

  !show_missing && return length(notdonecombs)

	for c in notdonecombs

		println(c)
		println()

	end
	
  return notdonecombs
	
end



function existing_data(task; show_existing=false)

	allcombs = task.get_paramcombs()

	donecombs = filter(c -> task.files_exist(c...), allcombs)

	println(string("\nTask: ",get_taskname(task),
								 "\nTotal number of jobs done: ",
									length(donecombs),"/",length(allcombs),"\n"))

  !show_existing && return length(donecombs)

	for c in donecombs

		println(c)
		println()

	end
	
end



#function delete_data(task; show_deleted=false)

#	allcombs = task.get_paramcombs()

#	Utils.Delete_NamesVals(


#end 



#===========================================================================#
#
#	Compute for one parameter combination (the first -- for testing purposes)
#
#---------------------------------------------------------------------------#

function get_data_one(task, pick::Function=Utils.DictFirstVals; 
											kwargs...)

	P = task.get_paramcombs(;repl=(l::Int,P...)->pick(P[l]))[1]

	return task.get_data(P...; force_comp=true, mute=true, kwargs...)

end



function get_plot_one(task, pick=Utils.DictFirstVals)

	good_P = task.get_paramcombs(;repl=(l::Int,P...)->pick(P[l]))[1]

	plot_P = task.get_plotparams(good_P...)

	return task.plot(plot_P)

end



#===========================================================================#
#
#	Compute for random parameters combinations (when there are too many combs)
#
#---------------------------------------------------------------------------#


function get_data_random(task; seed=nothing, check_data=true, mute=true)

	seed isa Int && Random.seed!(seed)

	for i=1:10^6

#		println(string(gethostname(),"\t",Utils.get_arg(1,2,Int64),"\t",i,"\n"))
		!mute && println(string(gethostname(),"\t",i,"\n"))
	
		Tasks.Generic.get_data_one(task, Utils.DictRandVals;
																			 force_comp=!check_data, mute=mute)
	
	end

end



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#









function choose_obs_i(data; P=nothing, k=nothing, f="last")

	getkeys,getvals,valofkey = Utils.gk_gv_vk_dictorJLDAW(data)


	if !isnothing(k) 

		v = valofkey(data, k)

		!isnothing(v) && return (k, v)

		error("Key '$k' not found. Try ",getkeys(data))

	end 


	K = sort(collect(getkeys(data)))

	
	i = if P isa Number 

				P

			elseif P isa AbstractString

				tryparse(Float64, P)
				
			elseif P isa NamedTuple

				get(P, :obs_i, nothing) 

			elseif !isnothing(P)

				Dict(P) |> function aux(p)

										haskey(p, "obs_i") && return p["obs_i"]

										return get(p, :obs_i, nothing)

									end 
			else 

				nothing

			end

	isnothing(i) || return choose_obs_i(data;k=get(K,Int(trunc(i)),nothing),f=f)


	f == "first" && return choose_obs_i(data; k=K[1])

	f == "last" && return choose_obs_i(data; k=K[end])
	
	f == "sum" && return (join(K,"+"), sum(getvals(data)))

	error("'f' must be 'last' or 'sum' or 'first'")

end



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function combine_files_exist(tasks::AbstractVector{CompTask})::Function

	function files_exist(args...; kwargs...)::Bool

		for t in tasks 
		
			t.files_exist(args...; kwargs...) || return false 

		end 

		return true 

	end

end 


function combine_get_data(tasks::AbstractVector{CompTask})::Function

	function get_data(P...; samplevectors=12345, kwargs...)::Vector

		if !(isa(samplevectors,Int) && samplevectors==12345)

					error("kwarg disabled")
		end 
#			return Any[SampleVectors(reshape(obs,:,1), P), I...]

		return [t.get_data(P...; kwargs...) for t in tasks]

	end  

end 


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#
#
#function restrict_number_subplots(d::Int, v::AbstractVector; max_nr_plotx=nothing, max_nr_ploty=nothing, kwargs...)
#
#	if d in [3,4] 
#		
#		max_nr = Dict(3=>max_nr_plotx, 4=>max_nr_ploty)[d] 
#
#		if !isnothing(max_nr) 
#		
#			N = Int(round(max_nr))
#	
#			length(v)>N && return v[Utils.RescaleInds(N, v, 1)]
#
#		end 
#
#	end 
#
#	return v
#
#end
#

function fillvals(dims::AbstractVector{<:AbstractVector{Int}},
									args...)::Vector{Any} 
		
	@assert length(args)==length(dims) 

	v = Vector{Any}(undef, maximum(maximum, dims))

	for (dim, arg) in zip(dims, args)
		
		test = Utils.isList(arg) && length(dim)==length(arg)  

		for (i,d) in enumerate(dim) 
			
			v[d] = test ? arg[i] : arg

		end

	end 

	return v

end 


function correct_ndims(expected_ndims::Int,										 
											 d::Union{Number, AbstractArray{<:Number}})

	expected_ndims==0 && return length(d)==1 ? d[1] : error("Wrong dim")

	nr_singletons = ndims(d) - expected_ndims 

	if nr_singletons==0

		return d 

	elseif nr_singletons > 0

		kept_dims = setdiff(1:ndims(d), findall(size(d).==1)[1:nr_singletons])

		return reshape(d, size(d)[kept_dims]...)

	else 

		@warn "Incorrect value" 

		return reshape(d, fill(1, abs(nr_singletons))..., size(d)...)

	end 

end 



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function getf_setval(dims::AbstractVector{<:AbstractVector{Int}})::Function

	dim = maximum(maximum, dims)

	dim<=2 && return function sv!(Z::AbstractArray, i::Tuple,
																A::Union{Number, AbstractArray{<:Number}}
																)::Nothing

		Z[fillvals(dims, i, :)] = A 

	end 


	function setval!(Z::AbstractArray, I::Tuple, a::Number)::Nothing 

		Z[I[1:dim-2]...][I[dim-1:dim]...] = a

		return 

	end 


	function setval!(Z::AbstractArray, i::Tuple,
									 A::AbstractArray{<:Number})::Nothing 

		I1, I2 = getindex.([fillvals(dims, i, axes(A))], [1:dim-2, dim-1:dim])

		for i1 in Base.product(I1...)
	
			Z[i1...][I2...] = A[get.([i1], dims[2], [:])...]
		
		end 

		return 

	end 

	return setval!

end  

#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#
#
#
#function getf_getinds(dims::AbstractVector{<:AbstractVector{Int}})::Function
#
#	dim = maximum(maximum, dims)
#
#	dim<=2 && return function gi(i::Tuple, 
#															 A::Union{Number, AbstractArray{<:Number}}
#															 )
#
##		fillvals(dims, i, axes(A))
#
#	end 
#
#
#	function gi(I::Tuple, a::Number)
#
#		I[1:dim-2] => I[dim-1:dim]
#
#	end 
#
#
#	function gi(i::Tuple, A::AbstractArray{<:Number})
#
#		I1, I2 = getindex.([fillvals(dims, i, axes(A))], [1:dim-2, dim-1:dim])
#
#		return [i1 => I2 for i1 in Base.product(I1...)]
#
#	end 
#
#	return gi 
#
#end  




#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#



function init_multitask(M, 
												internal_keys::AbstractVector, 
												ext_dim::AbstractVector=[]; 
												kwargs...) 
	#::AbstractVector{Pair{Int,<:AbstractVector{<:Number}}}=[]; kwargs...)


	rmv, add = Parameters.rmv_add_summarize(; kwargs...)

	rmv_ikey = Parameters.combine_functions_addrem(rmv,

							 map(collect(internal_keys)) do ik 

									Parameters.replace_parameter_fs(0, ik...)[1]

								end 
								)

	function add_ip(intpar::Utils.List)::Function 

		Parameters.combine_functions_addrem(add) do (level, P...)#::Int, P::Vararg)

			newp  = map(zip(internal_keys, intpar)) do ((k,ls),v)

				minimum(ls)==level ? k=>v : nothing

			end 

			return Utils.adapt_merge(P[l], filter(!isnothing, newp))

		end 

	end 



	if any(ext_dim.<2) && any(ext_dim.>=2) 
		
		@warn "Might have trouble constructing Z"

	end 


#	external_dims = [k for (k,v) in ext_par]

#	external_param = map(ext_par) do (k,v)
#
#		v
#		restrict_number_subplots(k, collect(v); kwargs...)
#		
#	end 



	return init_multitask_(M,
												 OrderedDict(internal_keys),
#												 [k for (k,v) in ext_par],
												ext_dim,												 
												 rmv_ikey,
												 add_ip,
													)


end 


	


function init_multitask_(M, 
												 internal_keys::OrderedDict,#{Symbol,Int},
												 external_dims::AbstractVector{Int},
												 rmv_internal_key::Function,
												 init_add_internal_param::Function,
												 )

	dim = length(internal_keys) + length(external_dims)

	internal_dims = setdiff(collect(1:dim), external_dims)

	dims = [internal_dims, external_dims]



	# -------------- names -------------- # 
	
	all_names = internal_names, external_names = [["x","y","plotx","ploty"][d] for d in dims]



	zname = filter(!in(vcat(all_names...)), ["x","y","z"])[1]

	internal_labels = [replace(string(k),"_"=>" ") for (k,v) in internal_keys]



	# ------- all internal params and their combinations ---- #

	internal_allparams = map(zip(internal_keys, internal_dims)) do ((k,l),d)

		v = Parameters.get_allP(M, fill(Dict(), minimum(l)-1)...)[k]

		return sort(collect(v))
#		return restrict_number_subplots(d, sort(v); kwargs...)

	end 


	combs(q) = collect(Base.product(q...))[:]
	
	internal_params_inds = combs(axes.(internal_allparams, 1))

					#-----------------#

	tasks = map(combs(internal_allparams)) do intpar  

						CompTask(M; 
										 rmv_internal_key = rmv_internal_key,
										 add_internal_param = init_add_internal_param(intpar),
														 )
					 end	
					 

	task0 = CompTask(M; rmv_internal_key=rmv_internal_key)

	@show dim 
@show dims
	@show all_names
	@show zname internal_labels internal_allparams internal_params_inds 


				# ----------------- #
				#
	setval! = getf_setval(dims)
#=

					#-----------------#


	function construct_Z(obs::AbstractString, P...; 
											 kwargs...)::Dict{String,Any}
		
		Data = get_data(P...; fromPlot=true, kwargs..., target=obs)

		d1 = Data[1][obs]

		if !Utils.is_dict_or_JLDAW(d1) 
			
			return merge!(construct_Z(D->D[obs], Data),
																			 Dict(zname*"label" => obs))
		end


		sub_obs = choose_obs_i(d1; P=P, f="first")[1]

		return merge!(construct_Z(D->choose_obs_i(D[obs]; k=sub_obs)[2], Data),
									Dict(zname*"label" => "$obs $sub_obs"),)

	end 


	
	function construct_Z(get_obs::Function, P...;
											 label=nothing,
											 kwargs...)::Dict{String,Any}
		
		construct_Z(get_obs, get_data(P...; fromPlot=true, kwargs...), label)

	end 


=#	 

#get_obs 

#external_param: info on what get_obs returns 
# length(external_param) == ndims(get_obs(...))
# dim==1 "x"  

#	dim==1 => return one vector 
#	dim==2 => return one matrix 
# dim==3 => return a vector of matrices 
# dim==4 => return a matrix of matrices 



#
#
#
#function setobs!(Z, inds, obs)
#	
#	if isempty(external_param) 
#	
#		Z[inds...] = obs[1]
#
#	else 
#
#		I = fillvals(dims, inds, :)
#		
#		Z[I...] = obs[axes.(external_param, 1)...] 
#
##		reshape(obs, :)
#
#	end 
#
#
#end


#function result(data)
#
##	dim in [1,2]
##
##	d1 = correct_ndims(data[1])
##
##	out = zeros(fillvals(dims, length.(internal_allparams), size(d1))...)
##
##	setval!(out, ?, d1)
##
##	---- 
#
#
#for I in allinds  
#
#
#I[1:dim-2]
#
#canvas, pixels = Utils.zipmap(zip(all_inds,data)) do (i,d)
#
#	s = fillvals(dims, length.(internal_allparams), size(d))
#
#	I = fillvals(dims, i, axes(d))
#
#
#	return (s[1:dim-2],s[dim-1:dim])
#
#end 
#
#
#@assert length(unique(canvas))==1 
#
#
#
#
##	n-2 
#
#println() 
#

#
#end 
#
#L = fillvals(dims, length.(internal_allparams),
#						 									length.(external_param))

#
#dim in [1,2] 


# d: one instance of the data 
#
# d = correct_ndims(d)  size-corrected data

#	s = fillvals(dims, length.(internal_allparams), size(d))) 
#

	println(fillvals(dims, length.(internal_allparams), 20))




	function construct_Z(get_obs::Function, Data::AbstractVector,
											 label::Nothing=nothing)::Dict{String,Any}
											 
		out = Dict{String,Any}(
							zname => zeros(fillvals(dims,
																			length.(internal_allparams),
																			length.(external_param))...))

		for (inds,D) in zip(internal_params_inds, Data)
			
			setobs!(out[zname], inds, get_obs(D))

		end 

		return out 

	end

#=
	function construct_Z(get_obs::Function, Data::AbstractVector,
											 label::AbstractString)::Dict{String,Any}

		out = construct_Z(get_obs, Data)
	
		if !isempty(label)

			out[zname*"label"] = label

		end 

		return out 

	end 





	# --------------- #

	nrowcol = map(3:dim) do i

	  pj = argmax([any(d.==i) for d in dims])

		return (internal_allparams, external_param)[pj][argmax(dims[pj].==i)]

	end  

=#

#construct_Z = nothing

	out_dict = merge(

#			Dict(zip(external_names,external_param)),
			
			Dict(zip(internal_names,internal_allparams)),

			Dict(zip(internal_names.*"label",internal_labels)),

									)




	return (CompTask(task0.name,
									 task0.get_plotparams, 
									 task0.get_paramcombs, 
									 combine_files_exist(tasks),
									 combine_get_data(tasks)
									 ),
					out_dict, 
					construct_Z,
#					nrowcol,
					)
end









#############################################################################
end

