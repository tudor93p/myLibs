module ComputeTasks
#############################################################################

using Distributed

import ..Random, PyCall, ..OrderedDict 

import ..Utils 
import ..Parameters





#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


get_taskname(t) = string(Base.fullname(t)[end])




#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

struct CompTask

	get_plotparams::Function 
	
	get_paramcombs::Function 
	
	files_exist::Function 
	
	get_data::Function 
	
end 




#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#




function init_task(M;

									found_files = M.FoundFiles,

								  read_data = M.Read,

					 			 	calc_data = M.Compute,

									valid_pcomb = nothing,
					
									kwargs...

									)


	rmv_internal_key, add_internal_param = Parameters.rmv_add_summarize(; kwargs...)


	fill_internal(P) = Parameters.fillParams_internalP(M, P, add_internal_param)




	files_exist_1(p,target) = found_files(p...; target=target)
	

	pr(d::Parameters.UODict) = println(Utils.NT(d))
	pr(d::Utils.List) = foreach(pr,d) 

	function calc_or_read(p, force_comp, mute, target)
	
		F = (!force_comp && files_exist_1(p, target)) ? read_data : calc_data

		mute && return F(p...; target=target)

		println("\n", Base.nameof(F), " ", Parameters.tostr(M) , "\n")
		
		pr(p)	
	
		out = @time F(p...; target=target)
	
	  println("Finished.\n")
	
	  return out
	
	end
	


	return CompTask(

		Parameters.f_get_plotparams(M, rmv_internal_key),
	
	#	function get_plotparams(P=nothing)
	
	#		Parameters.convertParams_toPlot(M, P; repl=rmv_internal_key)
	
	#	end,
	
	
		function get_paramcombs(;cond=nothing, repl=nothing)
	
			Parameters.get_paramcombs(M;
																
				repl = (rmv_internal_key, repl),
	
				cond = (valid_pcomb, cond),
				
				)
		end,
	
	
	
		function files_exist(P...; target=nothing)
	
			files_exist_1(fill_internal(P), target)
			
		end,
	
	
	
		function get_data(P...; target=nothing, force_comp=false,
											mute=true, shuffle=false, fromPlot=false,
											get_good_P=false, apply_rightaway=(d,p)->d)
	
			good_P = fill_internal(
									fromPlot ? Parameters.convertParams_fromPlot(M,P) : P
									)
	
			data = if target=="None" nothing else 
	
					apply_rightaway(calc_or_read(good_P, force_comp, mute, target),
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
														check_data=true, mute=true, )

	AC = task.get_paramcombs()

	check_data && filter!(!task.files_exist, AC)


	rev && reverse!(AC)


  if shuffle
	
		seed isa Int && Random.seed!(seed)

		Random.shuffle!(AC)

	end

  Utils.Distribute_Work(AC, task.get_data;
												force_comp=!check_data, mute=mute)#, shuffle=shuffle)

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

#	notdonecombs = filter(!task.files_exist, allcombs)

#	nprocs()>1 && Random.shuffle!(allcombs)

	id = show_missing ? isequal(first(workers())) : x->false

	notdonecombs = allcombs[
													
		(nprocs()>1 ? pmap : map)(enumerate(allcombs)) do (i,c)

			out = !task.files_exist(c)


			id(myid()) && println(i,"/",length(allcombs)," Files ", out ? "don't " : "", "exist")

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

	donecombs = filter(task.files_exist, allcombs)

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

function get_data_one(task, pick=Utils.DictFirstVals; kwargs...)

	task.get_data(task.get_paramcombs(;repl=(l,p)->pick(p))[1];
								force_comp=true, mute=true, kwargs...)

end



function get_plot_one(task, pick=Utils.DictFirstVals)

	p = task.get_paramcombs(;repl=(l,P)->pick(P))

	return p |> first |> task.get_plotparams |> task.plot 


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





function gk_gv_vk_dictorJLDAW(T)

	T<:AbstractDict && return keys,values,(d,k)->get(d,k,nothing)

	function valofkey(d,k)

		i = findfirst(isequal(k), d.keys)

		return !isnothing(i) ? d.values[i] : nothing 

	end 

	return (d->d.keys, d->d.values, valofkey) 

end 




function choose_obs_i(data::T; P=nothing, k=nothing, f="last") where T

	getkeys,getvals,valofkey = gk_gv_vk_dictorJLDAW(T)


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

				Dict(P) |> function (p)

										haskey(p, "obs_i") && return p["obs_i"]

										return get(p, :obs_i, nothing)

									end 
			else 

				nothing

			end

	if !isnothing(i) 
		
		return choose_obs_i(data; k=get(K, Int(trunc(i)), nothing), f=f)

	end 


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

function combine_files_exist(tasks::AbstractVector{CompTask})

	function files_exist(args...; kwargs...)

		for t in tasks 
		
			t.files_exist(args...; kwargs...) || return false 

		end 

		return true 

	end

end 


function combine_get_data(tasks::AbstractVector{CompTask})

	function get_data(P::AbstractDict; samplevectors=12345, kwargs...)

		if !(isa(samplevectors,Int) && samplevectors==12345)

					error("kwarg disabled")
		end 
#			return Any[SampleVectors(reshape(obs,:,1), P), I...]

		return [t.get_data(P; kwargs...) for t in tasks]

	end  

end 


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function restrict_number_subplots(d::Int, v::AbstractVector; max_nr_plotx=nothing, max_nr_ploty=nothing, kwargs...)

	max_nr = [max_nr_plotx, max_nr_ploty][d.==[3,4]]

	if !isempty(max_nr) && !isnothing(max_nr[1]) 
		
		N = Int(round(max_nr[1]))
	
		if length(v)>N
			
			inds = Utils.Rescale(1:N, axes(v,1)) .|> trunc .|> Int 

			return v[inds]

		end 

	end 

	return v

end


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#



function init_multitask(M, internal_keys_, ext_par=[]; kwargs...)

	rmv_internal_key, add_internal_param = Parameters.rmv_add_summarize(; kwargs...)

	internal_keys = OrderedDict(internal_keys_)

	external_dims = [k for (k,v) in ext_par]

	external_param = map(ext_par) do (k,v)

		restrict_number_subplots(k, collect(v); kwargs...)
		
	end 
	

	dim = length(internal_keys) + length(external_dims)

	internal_dims = filter(!in(external_dims),1:dim)

	dims = [internal_dims, external_dims]


	function fillvals(int, ext)

		v = Vector{Any}(undef, dim)

		for dim_par in zip(dims, (int,ext)), (i,vi) in zip(dim_par...)
				
			v[i] = vi

		end

		return v

	end


	names, external_names = [["x","y","plotx","ploty"][d] for d in dims]



	zname = filter(!in([names;external_names]),["x","y","z"])[1]


	labels = map(vcat(keys(internal_keys)...)) do k
																		replace(string(k),"_"=>" ") end


	allparams = map(zip(internal_keys,internal_dims)) do ((k,l),d)

		v = M.allparams(fill(nothing, minimum(l)-1)...)[k]

		return restrict_number_subplots(d, sort(v); kwargs...)

	end 


	internal_params, internal_params_inds = [
		collect(Base.product(q...))[:] for q in [allparams, axes.(allparams,1)]]


					#-----------------#



	rmv_ikey = Parameters.combine_functions_addrem(rmv_internal_key,

							 map(collect(internal_keys)) do ik 

									Parameters.replace_parameter_fs(0, ik...)[1]

								end 
								)


	ik_levels = [minimum(vcat(v...)) for v in values(internal_keys)]

	tasks = map(internal_params) do intpar  

						ip_all = collect(zip(keys(internal_keys), intpar))

						newp(l,p) = merge(p, Dict(ip_all[ik_levels.==l]))

						add=Parameters.combine_functions_addrem(newp, add_internal_param)

						return init_task(M;
														 rmv_internal_key = rmv_ikey,
														 add_internal_param = add,
														 )
					 end	
					 

	task0 = init_task(M; rmv_internal_key=rmv_ikey)


				# ----------------- #


	function setobs!(Z, inds, obs)
		
		if isempty(external_param) 
		
			Z[inds...] = obs[1]

		else 

			I = fillvals(inds, fill(:, length(external_dims)))
	
			Z[I...] = obs[axes.(external_param, 1)...] 

		end 


	end
					#-----------------#

	function construct_Z(obs::AbstractString, P::AbstractDict; kwargs...)
		
		Data = get_data(P; fromPlot=true, kwargs..., target=obs)

		d1 = Data[1][obs]

		if !Utils.is_dict_or_JLDAW(d1) 
			
			return merge!(construct_Z(D->D[obs], Data),
																			 Dict(zname*"label" => obs))
		end


		sub_obs = choose_obs_i(d1; P=P, f="first")[1]

		return merge!(construct_Z(D->choose_obs_i(D[obs]; k=sub_obs)[2], Data),
									Dict(zname*"label" => "$obs $sub_obs"),)

	end 


	
	function construct_Z(get_obs::Function, P::AbstractDict, label::T=nothing; 
											 kwargs...) where T<:Union{Nothing,AbstractString}
		
		construct_Z(get_obs, get_data(P; fromPlot=true, kwargs...), label)

	end 


	
	function construct_Z(get_obs::Function, Data::AbstractVector, 
											 label::T=nothing) where T<:Union{Nothing,AbstractString}

		out = Dict{AbstractString,Any}(
							zname => zeros(fillvals(length.(allparams),
																			length.(external_param))...))

		if isa(label,AbstractString) && !isempty(label)

			out[zname*"label"] = label

		end 


		for (inds,D) in zip(internal_params_inds, Data)
			
			setobs!(out[zname], inds, get_obs(D))

		end 

		return out 

	end
		






	# --------------- #

	out_dict = merge(

			Dict(zip(external_names,external_param)),
			
			Dict(zip(names,allparams)),

			Dict(zip(names.*"label",labels)),

									)



	nrowcol = map(3:dim) do i

	  pj = argmax([any(d.==i) for d in dims])

		return (allparams, external_param)[pj][argmax(dims[pj].==i)]

	end  


	return (CompTask(task0.get_plotparams, 
									 task0.get_paramcombs, 
									 combine_files_exist(tasks),
									 combine_get_data(tasks)),
					out_dict, 
					construct_Z,
					nrowcol,
					)

end









#############################################################################
end

