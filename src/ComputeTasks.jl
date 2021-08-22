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





function CompTask(C::Parameters.Calculation;

									valid_pcomb = nothing,
					
									kwargs...

									)::CompTask

	#found_files, read_data, calc_data take as arguments one or more 'UODict'. 
	# To be used as f(P::UODicts) or f(P...)

	rmv_internal_key, add_internal_param = Parameters.rmv_add_summarize(; kwargs...)


	get_plotparams = Parameters.f_get_plotparams(C, rmv_internal_key)

	fill_internal = Parameters.f_fill_internal(C, add_internal_param)





	function calc_or_read(p::Utils.List, force_comp::Bool, mute::Bool; 
												kwargs...)

		F,text = if !force_comp && C.FoundFiles(p...; kwargs...)
			
					(C.Read,"Reading files")

				else 

					(C.Compute,"Calculating")

				end  


		mute && return F(p...; kwargs...)

		println()

		@info string(text, " ", Parameters.tostr(C))
		
		println()
		
	

		for d in p 
			
			Utils.isList(d) ? println.(Utils.NT.(d)) : println(Utils.NT(d))

		end 

	
		out = @time F(p...; kwargs...)
	
	  println("Finished.\n")
	
	  return out
	
	end
	



	return CompTask(

		get_taskname(C), 

		get_plotparams, 
	
	
		function get_paramcombs(;cond=nothing, repl=nothing)
	
			Parameters.get_paramcombs(C;
																
				repl = (rmv_internal_key, repl),
	
				cond = (valid_pcomb, cond),
				
				)
		end,
	

	
		function files_exist(Ps...; target=nothing)

				C.FoundFiles(fill_internal(Ps)...; target=target) 
			
		end,
	
	
	
		function get_data(P...; target=nothing, 
											force_comp::Bool=false,
											mute::Bool=true, 
											shuffle::Bool=false, 
											fromPlot::Bool=false,
											get_good_P::Bool=false, 
											apply_rightaway::Function=(d,p)->d)

			good_P = fill_internal(
									fromPlot ? Parameters.convertParams_fromPlot(C, P...) : P
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



function get_plot_one(task, pick::Function=Utils.DictFirstVals)

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
	
		get_data_one(task, Utils.DictRandVals;
																			 force_comp=!check_data, mute=mute)
	
	end

end



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#




parse_obs_i(::Nothing, ::AbstractVector, f::Nothing=nothing)::Nothing = nothing 

function parse_obs_i(::Nothing, K::AbstractVector, f::AbstractString
										)::Union{Nothing,Int} 

	get(Dict("first"=>1, "last"=>length(K)), f, nothing)

end 

function parse_obs_i(i::Int, K::AbstractVector, args...)::Union{Nothing,Int}

	in(i,axes(K,1)) ? i : parse_obs_i(nothing, K, args...)

end 


function parse_obs_i(x::Float64, args...)::Union{Nothing,Int}
	
	parse_obs_i(Int(trunc(x)), args...)

end 
	


function parse_obs_i(s::AbstractString, args...)::Union{Nothing,Int} 
	
	parse_obs_i(tryparse(Float64,s), args...)

end 


function parse_obs_i(P::NamedTuple, args...)::Union{Nothing,Int} 
	
	parse_obs_i(get(P,:obs_i,nothing), args...)

end 


function parse_obs_i(P::AbstractDict, args...)::Union{Nothing,Int} 
	
	parse_obs_i(get(P, "obs_i", get(P, :obs_i, nothing)), args...)

end 

function parse_obs_i(P::Pair, args...)::Union{Nothing,Int}
	
	parse_obs_i(Dict(P), args...)

end 


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#




#function choose_obs_i(data::AbstractDict, k, v::Nothing; kwargs...)
#
#	error("Key '$k' not found. Try ",getkeys(data))
#
#end 

#function choose_obs_i(data, k, v; kwargs...)::Tuple{Any,Any}
#
#	(v,k)
#
#end 

function choose_obs_i(data::AbstractDict{Tk,Tv}, k::Tk; 
											kwargs...)::Tuple{Tv, Tk} where {Tk,Tv}

	(data[k], k)
#	choose_obs_i(data, k, Utils.valofkey_dictorJLDAW(data, k))

end 



function choose_obs_i(data::AbstractDict{Tk,Tv},
											k::Nothing=nothing; 
											P=nothing, f=nothing, kwargs...
											)::Tuple{Tv,<:Any} where {Tk,Tv}


#	getkeys,getvals,valofkey = Utils.gk_gv_vk_dictorJLDAW(data)

#	K,i = sort(collect(getkeys(data)),by=string), parse_obs_i(P)
	K = sort(collect(keys(data)), by=string)
	i = parse_obs_i(P, K, f) 



	#	!isnothing(i) && in(i,axes(K,1)) && return choose_obs_i(data, K[i])

	isnothing(i) || return choose_obs_i(data, K[i])
	
	isa(f,AbstractString) && f=="sum" && return sum(values(data)), join(K,"+")
		
	error("Provide a valid 'k' or 'f'. 'f' must be 'last' or 'sum' or 'first'")

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

	v = Vector{Any}(undef, sum(length, dims))

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


function getf_setval(dims::AbstractVector{<:AbstractVector{Int}}
										 )::Tuple{Function,Function}


	dim = sum(length, dims)

	if dim<=2 
		
		function gi1(i::NTuple{N, Int} where N,
								A::Union{Number, AbstractArray{<:Number}}
								)::Vector{Union{Int,AbstractVector{Int}}}

			fillvals(dims, i, axes(A))

		end 
	

		function sv1!(Z::AbstractArray, i::Tuple{Vararg{Int}}, 
								A::Union{Number, AbstractArray{<:Number}}
									)::Nothing 

			Z[gi1(i, A)...] = A 

			return 

		end 


		return gi1,sv1! 

	end 




	function gi2(I::NTuple{N, Int} where N, a::Number
							)::Tuple{Tuple{Tuple{Vararg{Int}}}, 
											 Tuple{Vararg{Int}}}

		(I[1:dim-2],), I[dim-1:dim]

	end 


	function sv2!(Z::AbstractArray, I::NTuple{N, Int} where N, a::Number
									 )::Nothing 

		i1,i2 = gi2(I, a)

		Z[i1[1]...][i2...] = a

		return 

	end 




	function gi2(i::NTuple{N, Int} where N, A::AbstractArray{<:Number}
							)::Tuple{Base.Iterators.ProductIterator,
											 Vector{Union{Int, AbstractVector{Int}}}}
		
		I1, I2 = getindex.([fillvals(dims, i, axes(A))], [1:dim-2, dim-1:dim])

		return Base.product(I1...), I2

	end 


	function sv2!(Z::AbstractArray, i::NTuple{N, Int} where N,
									 A::AbstractArray{<:Number})::Nothing 

		I1, i2 = gi2(i, A)

		for i1 in I1 

#			Z[i1...][i2...] = A[get.([i1], dims[2], [:])...] 

			Z[i1...][i2...] = A[(ed<dim-1 ? i1[ed] : Colon() for ed in dims[2])...]
		
		end 

		return 

	end 

	return gi2,sv2!

end  

#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function initialize_dataZ(dims::AbstractVector{<:AbstractVector{Int}},
													inds::AbstractVector{<:Tuple{Vararg{Int}}}, 
													data::AbstractVector
													)::Array

	geti = getf_setval(dims)[1]

	dim = sum(length, dims)

	function update_max!(S::AbstractVector{Int}, 
											 vals::Union{AbstractVector{Int}, Tuple{Vararg{Int}}})

#		S .= max.(S,vals) 
#		return 

		for (i, (s,v)) in enumerate(zip(S, vals))

			v>s && setindex!(S, v, i)

		end  

	end 


	if dim<=2 

		shape = zeros(Int, dim) 
	
		for iA in zip(inds, data)

			update_max!(shape, maximum.(geti(iA...)))

		end 
	
		return zeros(shape...)

	end 


	shape = Dict{Tuple{Vararg{Int}}, Vector{Int}}()


	for iA in zip(inds, data)

		I1,I2 = geti(iA...) 

		i2 = [maximum(i2_) for i2_ in I2]

		for i1 in I1 

			if haskey(shape, i1) 

				update_max!(shape[i1], i2) 

			else 

				shape[i1] = i2

			end 

		end 

	end 

	# for each panel: values for external params within the panel 
	# external_dims...?
	#

	panels = zeros(Int, dim-2)

	for panel in keys(shape) 

		update_max!(panels, panel)

	end  




	return [zeros(shape[ci.I]...) for ci in CartesianIndices(Tuple(panels))]
	
end 


function fill_dataZ!(Z::AbstractArray,
										 dims::AbstractVector{<:AbstractVector{Int}},
													inds::AbstractVector{<:Tuple{Vararg{Int}}}, 
													data::AbstractVector
													)::Array

	setval! = getf_setval(dims)[2]

	for item in zip(inds, data) 

		setval!(Z, item...)
		
	end 

	return Z 

end 


function fill_dataZ(dims::AbstractVector{<:AbstractVector{Int}},
													inds::AbstractVector{<:Tuple{Vararg{Int}}}, 
													data::AbstractVector
													)::Array

	fill_dataZ!(initialize_dataZ(dims, inds, data), dims, inds, data)

end 

 



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


#dictJLDAW_getval(k)::Function = D -> Utils.valofkey_dictorJLDAW(D, k) 

function getval(k)::Function

	function getval_(D::AbstractDict{TK,TV})::TV where {TK,TV}
		
		D[k]

	end 

end 


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function parse_external_params(int_dim::Int,
		 ext_pairs_::AbstractVector{<:Pair{Int,<:Any}}=Pair{Int,Int}[],
			ext_lab_::AbstractVector{<:AbstractString}=fill("",length(ext_pairs_)),
												)::Tuple{Vector{Int}, 
																 Vector{Union{Function,
																							<:AbstractVector{<:Number}}},
																 Vector{String}}

	isempty(ext_pairs_) && return (Int[], Function[], String[])


	@assert length(ext_pairs_)==length(ext_lab_) 


	dim = int_dim + length(ext_pairs_)

	@assert 1<=dim<=4 # cannot plot more 

	order = sortperm(ext_pairs_, by=first)

	ext_dim, ext_par = collect.(zip(ext_pairs_[order]...))



	@assert allunique(ext_dim) & issubset(ext_dim, 1:dim)


	Ts = vs,ss,fs,Is = [isa.(ext_par,T) for T in [AbstractVector{<:Number},
																								AbstractString, 
																								Function, 
																								Int]]

	@assert sum(count, Ts)==length(ext_par) "Unsupported type"  

	dim>2 && @assert all(vs[ext_dim.<dim-1]) "Panel labels must be explicit."



	make_function(obs::AbstractString)::Function = getval(obs)
			

	function make_function(i::Int)::Function

		function get_plotaxis_(D::AbstractArray{<:Number})::Vector{Int}
			
			axes(D,i)

		end 
		
	end 



	return (ext_dim,
					[si ? make_function(v) : v for (v,si) in zip(ext_par, ss.|Is)],
					collect(ext_lab_[order]),
					)

																											




#
#	make_function(f::Function)::Function = f
#
#
#
#	function make_function(ext_par_method::AbstractVector{<:Real})
#
#		D -> ext_par_method
#
#	end 
#
#
#
#	any(s_or_f) && return make_function.(getindex.(ext_par, 2))
#
#	return getindex.(ext_par, 2)
#
#


end 


function assess_ext_par(Z::AbstractArray,
												dims::AbstractVector{<:AbstractVector{Int}},
												inds::AbstractVector{<:Tuple{Vararg{Int}}},
												Data::AbstractVector,
												zname::AbstractString,
												Keys::AbstractVector{<:AbstractString},
												candidates::AbstractVector{Union{Function,
																							<:AbstractVector{<:Number}}},
												zlab::Union{Nothing,AbstractString}=nothing,
												)::Dict{String,Any} 

	out = Dict{String,Any}()

	!isnothing(zlab) && !isempty(zlab) && setindex!(out, zlab, zname*"label")


	geti = getf_setval(dims)[1]


	for (K, f) in zip(Keys, candidates) 

		(f isa Function && applicable(f, Data[1])) || continue 
		
		if Z isa AbstractArray{<:Number} 
			
			out[K] = collect(f(Data[1]))

		elseif Z isa AbstractArray{<:AbstractArray{<:Number}}

			out[K] = similar(Z, Vector{Real})

			for (i, data) in zip(inds, Data)
				
				v = collect(f(data))

				for panel in unique(geti(i, data)[1])

					out[K][panel...] = v 

				end 
				
			end 

		end 

	end 

	return out 

end
	

function get_Z_and_axes(
												dims::AbstractVector{<:AbstractVector{Int}},
												inds::AbstractVector{<:Tuple{Vararg{Int}}},
												Data::AbstractVector,
												zname::AbstractString,
												ext_par_args...
												)::Dict{String,Any} 

	Z = initialize_dataZ(dims, inds, Data) 

	return Dict{String,Any}(zname=>fill_dataZ!(Z, dims, inds, Data),
													assess_ext_par(Z, dims, inds, Data, zname,
																				 ext_par_args...)...)
end 



function get_Z_and_axes(obs::AbstractString, args...)::Dict{String,Any}

	get_Z_and_axes(getval(obs), args...)

end 



function get_Z_and_axes(get_obs::Function, 
												dims::AbstractVector{<:AbstractVector{Int}},
												inds::AbstractVector{<:Tuple{Vararg{Int}}},
												Data::AbstractVector,
												zname::AbstractString,
												ext_par_args...
												)::Dict{String,Any} 

		NData = get_obs.(Data) # numerical data

		Z = initialize_dataZ(dims, inds, NData) 

		return Dict{String,Any}(zname=>fill_dataZ!(Z, dims, inds, NData),
														assess_ext_par(Z, dims, inds, Data, zname, 
																					 ext_par_args...)...,
#														assess_ext_par(Z, dims, inds, NData, zname, 
#																					 ext_par_args...)...,
														)
end 


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#




function init_multitask(C::Parameters.Calculation,
												internal_keys::AbstractVector, 
												ext_args...;
												kwargs...) 
	#::AbstractVector{Pair{Int,<:AbstractVector{<:Number}}}=[]; kwargs...)

	rmv, add = Parameters.rmv_add_summarize(; kwargs...)

	rmv_ikey = Parameters.combine_functions_addrem(rmv,

			 			map(collect(internal_keys)) do (key,level)

									Parameters.replace_parameter_fs(nothing, key, level)[1]

								end 
								)


	function add_ip(intpar::Utils.List)::Function 

		function add_(level::Int, P...)
	
			inds = [minimum(ls)==level for (k,ls) in internal_keys]

			any(inds) || return P[level]

			return Utils.adapt_merge(P[level], zip(first.(internal_keys[inds]),
																						 intpar[inds]))

		end 
		
		return Parameters.combine_functions_addrem(add, add_)

	end 



#
#
#	intpar[inds]
#
#	replace_parameter_fs((args...)->intpar[inds],
#											 first.(internal_keys[inds]),
#											 level)[2](l, P...)




#	external_param = map(ext_par) do (k,v)
#
#		v
##		restrict_number_subplots(k, collect(v); kwargs...)
#		
#	end 
#


	return init_multitask_(C,
												 OrderedDict(internal_keys),
												 parse_external_params(length(internal_keys),
																							 ext_args...)...,
												 rmv_ikey,
												 add_ip,
													)


end 


	


function init_multitask_(C::Parameters.Calculation, 
												 internal_keys::OrderedDict,#{Symbol,Int},
												 external_dims::AbstractVector{Int},
												 get_plotaxes::AbstractVector{<:Union{Function, <:AbstractVector{<:Number}}},
												 external_labels::AbstractVector{<:AbstractString},
												 rmv_internal_key::Function,
												 init_add_internal_param::Function,
												 )

	dim = length(internal_keys) + length(external_dims)

	internal_dims = setdiff(collect(1:dim), external_dims)

	dims = [internal_dims, external_dims]



	# -------------- names -------------- # 



#	dim_order = ["x","y","plotx","ploty"] # obsolete !

	dim_order = [["plotx","ploty"][1:dim-2];["x","y"][1:min(dim,end)]]

	all_names = internal_names, external_names = [dim_order[d] for d in dims]


	zname = filter(!in(vcat(all_names...)), ["x","y","z"])[1]





	# ------- all internal params and their combinations ---- #

	internal_allparams = map(zip(internal_keys, internal_dims)) do ((k,l),d)

		v = Parameters.get_allP(C, fill(Dict(), minimum(l)-1)...)[k]

		return sort(collect(v))
#		return restrict_number_subplots(d, sort(v); kwargs...)

	end 


	combs(q) = collect(Base.product(q...))[:]
	
	internal_params_inds = combs(axes.(internal_allparams, 1))

					#-----------------#

	tasks = map(combs(internal_allparams)) do intpar  

						CompTask(C; 
										 rmv_internal_key = rmv_internal_key,
										 add_internal_param = init_add_internal_param(intpar),
														 )
					 end	
					 


					 # get_plotparams and get_paramcombs will not always work; some keys removed at lower levels might be needed at higher levels 

	multitask = CompTask(getproperty.([CompTask(C; rmv_internal_key=rmv_internal_key)], [:name, :get_plotparams, :get_paramcombs])...,
									 combine_files_exist(tasks),
									 combine_get_data(tasks)
									 )

				# ----------------- #


#	P = rand(task0.get_paramcombs())

#	println() 

#	data = [t.get_data(P...) for t in tasks]

#	Z = initialize_dataZ(dims, internal_params_inds, data) 

#	Z = fill_dataZ(dims, internal_params_inds, data)

#=


#get_obs 

#external_param: info on what get_obs returns 
# length(external_param) == ndims(get_obs(...))
# dim==1 "x"   


#	dim==1 => return one vector 
#	dim==2 => return one matrix 
# dim==3 => return a vector of matrices 
# dim==4 => return a matrix of matrices 
					#-----------------#



=#	 

#	function zname_label(label::Nothing=nothing)::Dict{String,Any}
#		
#		Dict{String,Any}()
#
#	end
#
#	function zname_label(label::AbstractString)::Dict{String, Any}
#		
#		isempty(label) && return zname_label()
#
#		return Dict{String,Any}(zname*"label"=>label)
#
#	end 


	# --------- arrange data into the desired structure --------- #


	function construct_Z(Data::AbstractVector{<:Union{Number, <:AbstractArray{<:Number}}},
											 label...)::Dict{String,Any}

		get_Z_and_axes(dims, 
									 internal_params_inds, 
									 Data, 
									 zname, 
									 external_names, 
									 get_plotaxes, 
									 label...)


	end 


	function construct_Z(obs::Union{AbstractString, Function}, 
											 Data::AbstractVector, label...)::Dict{String,Any}
		
		get_Z_and_axes(obs,
									 dims, 
									 internal_params_inds, 
									 Data, 
									 zname, 
									 external_names, 
									 get_plotaxes, 
									 label...)
													
	end 





	# ----- called by the plot functions ---- # 


	function construct_Z(P::UODict, label...; kwargs...)::Dict{String,Any}

		construct_Z(multitask.get_data(P; fromPlot=true, kwargs...), 
								label...)

	end 


	function construct_Z(get_obs::Function, P::UODict, label...;
											 kwargs...)::Dict{String,Any}
		
		construct_Z(get_obs, 
								multitask.get_data(P; fromPlot=true, kwargs...), 
								label...)

	end 


	function construct_Z(obs::AbstractString, P::UODict; 
											 kwargs...)::Dict{String,Any}
		
		Data = multitask.get_data(P; fromPlot=true, kwargs..., target=obs)


		d1 = Data[1][obs]


		# if it's an array or number, not dictionary-like object

#		Utils.is_dict_or_JLDAW(d1) || return construct_Z(obs, Data)
		isa(d1, AbstractDict) || return construct_Z(obs, Data)
	


		sub_obs = choose_obs_i(d1; P=P, kwargs...)[2] 

		return construct_Z(D->D[obs][sub_obs], Data, "$obs $sub_obs") 

	end 


	# --------------- #
 
	out_dict = Dict{String,Any}(zip(internal_names,internal_allparams))


	for (n,(k,v)) in zip(internal_names, internal_keys)
	
		out_dict[n*"label"] = replace(string(k),"_"=>" ")
	
	end 
		

	for (n,a,k) in zip(external_names, get_plotaxes, external_labels)
		
		a isa AbstractVector{<:Number} && setindex!(out_dict, collect(a), n)

		out_dict[n*"label"] = replace(string(k), "_"=>" ") 

	end 


	return (multitask,
					out_dict, 
					construct_Z,
					[length(out_dict[k]) for k in ["plotx","ploty"][1:dim-2]]
					)
end









#############################################################################
end

