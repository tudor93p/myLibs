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



		mute && return F(p...; kwargs..., check_data=!force_comp)

		println()

		@info string(text,": ", Parameters.tostr(C))
		
		println()
		
	

		for d in p 
			
			Utils.isList(d) ? println.(Utils.NT.(d)) : println(Utils.NT(d))

		end 

	
		out = @time F(p...; kwargs..., check_data=!force_comp)

	  println("Finished.\n")
	
	  return out
	
	end



	function prep_out(data, good_P::AbstractVector, get_good_P::Bool, f::Nothing=nothing)

		get_good_P ? (data, good_P) : data 

	end 


	function prep_out(data, good_P::AbstractVector, 
										get_good_P::Bool, f::Function) 
		
		prep_out(f(data, good_P), good_P, get_good_P)

	end 





	function get_data(P::Union{Tuple,AbstractVector}; 
										fromPlot::Bool=false,
										target=nothing, 
										force_comp::Bool=false,
										mute::Bool=true, 
										get_good_P::Bool=false, 
										apply_rightaway=nothing,
										kwargs...)

		@assert !fromPlot  


		good_P = fill_internal(P)  

		target=="None" && return prep_out(nothing, good_P, get_good_P)

		data = calc_or_read(good_P, force_comp, mute; target=target)

		return prep_out(data, good_P, get_good_P, apply_rightaway)

	end 


	function get_data(P::AbstractDict; fromPlot::Bool=false, kwargs...)

		fromPlot || return get_data([P]; fromPlot=false, kwargs...)

		return get_data(Parameters.convertParams_fromPlot(C, P); 
										fromPlot=false, kwargs...) 

	end  

	function get_data(P...; fromPlot::Bool=false, kwargs...)

		@assert !fromPlot 

		return get_data(P; fromPlot=false, kwargs...)

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
	
	
		get_data, 	


	

 )

end

					






#===========================================================================#
#
#	Compute for all parameter combinations 
#							(optional: for which there is no data already stored)
#
#---------------------------------------------------------------------------#


function get_data_all(task; 
											shuffle::Bool=false, seed=nothing, 
											rev::Bool=false,
											check_data::Bool=true,
											mute::Bool=true, kwargs...)


	AC = task.get_paramcombs()

	check_data && filter!(c->!task.files_exist(c...; kwargs...), AC)

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

function missing_data(task; show_missing::Bool=false, kwargs...)

	allcombs = task.get_paramcombs()

	id = show_missing ? isequal(first(workers())) : x->false

	notdonecombs = allcombs[
													
		(nprocs()>1 ? pmap : map)(enumerate(allcombs)) do (i,c)


			out = !task.files_exist(c...; kwargs...)


#			id(myid()) && println(i,"/",length(allcombs)," Files ", out ? "don't " : "", "exist")

			return out

		end ]




	println(string("\nTask: ",get_taskname(task),
								 "\nTotal number of jobs left to do: ",
									length(notdonecombs),"/",length(allcombs),
									"=", round(100*length(notdonecombs)/length(allcombs),
														 digits=1),
									"%\n"))

  !show_missing && return length(notdonecombs)

	for c in notdonecombs

		println(c,"\n")

	end
	
  return notdonecombs
	
end



function existing_data(task; show_existing::Bool=false, kwargs...)

	allcombs = task.get_paramcombs()

	donecombs = filter(c -> task.files_exist(c...; kwargs...), allcombs)

	println(string("\nTask: ",get_taskname(task),
								 "\nTotal number of jobs done: ",
									length(donecombs),"/",length(allcombs),
									"=", round(100*length(donecombs)/length(allcombs),
														 digits=1),
									"%\n"))

  !show_existing && return length(donecombs)

	for c in donecombs

		println(c,"\n")

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

function choose_obs_i(data::T, args...; kwargs...
											)::Tuple{T, Nothing} where T<:Union{<:Number,<:AbstractArray{<:Number}}

	(data, nothing)

end 

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

	if isa(f,AbstractString) 
		
		f=="sum" && return sum(values(data)), join(K,"+")

		if f=="diff"

			@assert length(K)==2 

			k1,k2 = K 

			return data[k1]-data[k2], "$k1-$k2"

		end 

	end 
		
	error("Provide a valid 'k' or 'f'. 'f' must be 'last' or 'sum' or 'diff' or 'first'")

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


function combine_get_data(tasks::AbstractVector{Tuple{CompTask,String}},
#													NTuple{N,CompTask},
#													print_intpars::NTuple{N, AbstractString},
#													internal_keys::AbstractVector{Symbol},
#													internal_params#::OrderedDict=OrderedDict()
													)::Function #where N

	max_len = maximum([length(pr_ip) for (t,pr_ip) in tasks]) 
#	max_len = maximum(length, print_intpars) 

	max_len += 5

	function get_data(P...; samplevectors=12345, kwargs...
										)::Vector

		@assert isa(samplevectors,Int)&&samplevectors==12345 "kwarg disabled"

		pr = !haskey(kwargs, :mute)  

		pr && println()

		out = map(tasks) do (t,s)

			pr && print("\r$s",repeat(" ",max_len-length(s))) 

			return t.get_data(P...; kwargs...)

		end  

		pr && println("\r",repeat(" ",max_len))
		
		return out 

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
									args...)::Vector
		
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
								A::AbstractArray{<:Number}
								)::Vector#{<:Union{Int,AbstractVector{Int}}}

			out = fillvals(dims, i, axes(A)) 

			for out_i in out 

				@assert isa(out_i,Int) || isa(out_i,AbstractVector{Int})

			end 

			return out
			#{<:Union{Int,AbstractVector{Int}}}

		end 
	
		function gi1(i::NTuple{N, Int} where N,
								A::Number, 
								)::Vector

			gi1(i, [A])

		end 


		function sv1!(Z::AbstractArray, i::Tuple{Vararg{Int}}, 
								A::Union{Number, AbstractArray{<:Number}}
									)::Nothing 

			I = gi1(i, A) 

			sZ = size(view(Z, I...))

			sA = size(A) 

			if sZ==sA 

				Z[I...] = A 

			elseif prod(sZ)==prod(sA)

				Z[I...] .= A 

			else 

				error("Size mismatch: $sZ != $sA")

			end 

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
	

function get_Z_and_axes(dims::AbstractVector{<:AbstractVector{Int}},
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

	@assert allunique(internal_keys)

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

function add_line!(data::AbstractDict,
									 P::UODict,
									 k::Union{Symbol,AbstractString},
									 n::Char
									 )::Nothing 

	add_line!(data, P, k, string(n))

end 
function add_line!(data::AbstractDict,
									 P::UODict,
									 k::Union{Symbol,AbstractString},
									 n::Int
									 )::Nothing 

	add_line!(data, P, k, 'x'+n-1)

end 

function add_line!(data::AbstractDict,
									 P::UODict,
									 k::Symbol,
									 n::AbstractString
									 )::Nothing 

	in(n,["x","y"]) || return 
	
	for q in (k,string(k))

		q in keys(P) || continue 
		
		setindex!(data, P[q], n*"line")

		return 

	end 

	return 

end 
	
function add_line!(data::AbstractDict,
									 P::UODict,
									 k::AbstractString,
									 n::AbstractString
									 )::Nothing

	in(n,["x","y"]) || return 
	
	for q in (k,Symbol(k))

		q in keys(P) || continue 
		
		setindex!(data, P[q], n*"line")

		return 

	end 

	return 

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

	emptyP = [NamedTuple() for i=1:maximum(maximum, values(internal_keys))]


	tasks = map(combs(internal_allparams)) do intpar  

						add_ip = init_add_internal_param(intpar) 

						return (

						CompTask(C; 
										 rmv_internal_key = rmv_internal_key,
										 add_internal_param = add_ip,
														 ),

						string(merge((add_ip(i, emptyP...) for i=1:length(emptyP))...)),

						)

	end	
					 

	not_possib(args...; kwargs...) = error("The functions 'get_plotparams(P...)' and 'get_paramcombs()' do not always work for multitasks because some keys removed at lower levels might be needed by 'allparams' at higher levels. The two functions are not provided for multitasks.")

	all_plot_params = Parameters.convertParams_toPlot(C; repl=rmv_internal_key) 

	get_plotparams() = all_plot_params

	get_plotparams(arg1, args...) = not_possib()

	 
	multitask = CompTask(C.name, get_plotparams, 
											 not_possib, 
											 combine_files_exist(first.(tasks)),
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


	# --------- # 

	function data_and_startdict(P::UODict; kwargs...)

		sd = Dict{String,Any}()
		
		for (n,k) in zip(internal_names, string.(keys(internal_keys)))

			add_line!(sd, P, k, n)

			#	in(n,["x","y"]) && haskey(P, k) && setindex!(sd, P[k], n*"line")
	
		end 

		return multitask.get_data(P; fromPlot=true, kwargs...), sd

	end 




	function getval1(obs::AbstractString)::Function 

		getval_(data::AbstractDict) = data[obs]

		getval_((data, good_P)::Tuple{AbstractDict,AbstractVector}) = (data[obs], good_P)

		return getval_ 

	end 


	function getval2((data, good_P)::Tuple{T,AbstractVector}
									 )::T where T<:Union{AbstractDict, AbstractArray, Number}
		data 

	end 

	function getval2(data::T
									 )::T where T<:Union{AbstractDict, AbstractArray, Number}
		data 

	end  


	function extract_obs(d, obs::AbstractString; kwargs...)::Tuple{Function,String}

		extract_obs_(getval1(obs)(getval2(d)), obs; kwargs...)

	end 

	function extract_obs_(d::Union{Number, <:AbstractArray{<:Number}}, 
											 obs::AbstractString; kwargs...
											 )::Tuple{Function,String}
	
		getval1(obs), obs
	
	end 
	
	
	function extract_obs_(d::AbstractDict, obs::AbstractString; kwargs...
											 )::Tuple{Function,String}
	
		sub_obs = choose_obs_i(d; kwargs...)[2] 
	
		return getval1(sub_obs) ∘ getval1(obs), "$obs/$sub_obs"
	
	end 
	
	

	# ----- called by the plot functions ---- # 


	function construct_Z_(startdict::AbstractDict, args...)::Dict{String,Any}

		merge!(construct_Z(args...), startdict) ######### 

	end 



	function construct_Z(P::UODict, label...; kwargs...)::Dict{String,Any}

		Data, startdict = data_and_startdict(P; kwargs...)

		return construct_Z_(startdict, Data, label...) #########

	end 


	function construct_Z(get_obs::Function, P::UODict, label...;
											 kwargs...)::Dict{String,Any}
		
		Data, startdict = data_and_startdict(P; kwargs...)

		return construct_Z_(startdict, get_obs, Data, label...)

	end 




	function construct_Z(get_obs::Function, obs::AbstractString, P::UODict; 
											 kwargs...)::Dict{String,Any}
		
		Data, startdict = data_and_startdict(P; kwargs..., target=obs)

		f,lab = extract_obs(Data[1], obs; P=P, kwargs...)

		return construct_Z_(startdict, get_obs∘f, Data, lab)

	end 

	function construct_Z(obs::AbstractString, P::UODict; 
											 kwargs...)::Dict{String,Any}
	
		construct_Z(identity, obs, P; kwargs...) 

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

