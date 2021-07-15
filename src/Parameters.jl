module Parameters
#############################################################################

import ...OrderedDict, ...Utils 
import Random 

const ODict = Union{<:OrderedDict, <:NamedTuple}
const UDict = AbstractDict
const UODict = Union{<:AbstractDict, <:OrderedDict, <:NamedTuple}
const UODicts = Vararg{<:UODict}


#===========================================================================#
#
# typical allparams function
#
#---------------------------------------------------------------------------#

#function typical_allparams(M::Module, args...)
#
#	typical_allparams(M.usedkeys, args...)
#
#end

function typical_allparams(usedkeys::Function, args...)::Dict{Symbol,Any}

	typical_allparams(usedkeys(), args...)

end



function typical_allparams(usedkeys::AbstractVector{<:Symbol},
													 allp::UODict)::Dict{Symbol,Any}

	Dict(Symbol(k)=>allp[k] for k in intersect(keys(allp),usedkeys))

end








#===========================================================================#
#
# applies several add or remove functions recursively
#
#---------------------------------------------------------------------------#


function combine_functions_addrem(args...)::Function

	good_fs = Utils.flat(args...; keep=x->isa(x,Function))

	function combined_addrem(l::Int, p::UODict)

		p1 = copy(p)

		for fi in good_fs

			p1 = fi(l, p1)

		end 

		return p1

	end

	function combined_addrem(l::Int, p::Utils.List)::Vector
		
		combined_addrem.(l,p)

	end 


	return combined_addrem

end 




#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#



function combine_functions_cond(args...)::Function

	good_fs = Utils.flat(args...; keep=x->isa(x,Function))

	function combined_cond(P::UODict, l::Int)::Bool

		for cond in good_fs
			
			cond(P,l) || return false 

		end 

		return true 

	end


	function combined_cond(P::Utils.List, l::Int)::Bool 

		all(combined_cond.(P,l))

	end 


	return combined_cond


end 









#===========================================================================#
#
# Transforms the symbols into strings for python compatibility
# 	Since the system could contain components for which the same parameters
# 		 are specified (i.e. many leads), the label of each component is
# 		 incorporated in the new string key.
#
#---------------------------------------------------------------------------#

function convertKey_toPlot(p::UODict, k; separator="__")::String 

	haskey(p,:Label) && return join([p[:Label], k], separator) 

	return string(k)

end


function convertPair_toPlot(p::UODict, k; kwargs...)::Pair{String,<:Any}

	convertKey_toPlot(p, k; kwargs...) => p[k]

end 


function convertPairs_toPlot(p::UODict; kwargs...)::Vector{Pair{String,<:Any}}
	
	[convertPair_toPlot(p, k) for k in keys(p)] 

end 



##function merge_allparams(Ms::Vararg{Module,N}; args=[], kwargs...) where N
#function merge_allparams(inputs...; args=[], kwargs...)
#
#
#	out = map(inputs) do i
#	
#						typeof(i)<:Module && return i.allparams(args...; kwargs...)
#				
#						typeof(i)<:AbstractDict && return i
#																# already the result of typical_allparams
#			
#						return typical_allparams(i...)
#				
#					end
#					
##	Utils.isList(out,AbstractDict)
#
#	eltype(out) <: Union{NamedTuple,AbstractDict} && return merge(out...)
#
#	if Utils.isList(eltype(out),Union{NamedTuple,AbstractDict})
#
#		return [merge(o...) for o in zip(out...)]
#
#	end
#
#
#
#	error("Element type not undestood")
#
#end
#
#
#
#

#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function replace_parameter_fs(getP, Keys, level::Int=1)::Tuple{Function,Function}

	replace_parameter_fs(getP, Keys, [level])

end 

function replace_parameter_fs(getP::Function, Keys, 
															levels::AbstractVector{Int}
														 )::Tuple{Function,Function}

	K = Utils.flat(Keys) 

	rmv,aux = replace_parameter_fs(nothing, Keys, levels)


	function add(l::Int,p::T)::T where T<:UODict 

		in(l, levels) || return p
		
		return merge(p, Dict(zip(K, getP(l,p))))

	end 


	return rmv,add

end 



function replace_parameter_fs(getP, Keys, 
															levels::AbstractVector{Int}
														 )::Tuple{Function,Function}

	K = Utils.flat(Keys)

	function rmv(l::Int, p::T)::T where T<:UODict 

		in(l,levels) || return p 
		
		return filter(kv->!in(kv.first, K), p)

	end 


	return rmv,(l,p)->p

end 


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function rmv_add_summarize(; rmv_internal_key = nothing, 
													 	 add_internal_param = nothing,
														 constrained_params = nothing,
														 kwargs...)


	isnothing(constrained_params) && return (
								combine_functions_addrem(rmv_internal_key),
								combine_functions_addrem(add_internal_param)
																					)


	new_rmvs,new_adds = Utils.zipmap(constrained_params) do (level,cp)

			replace_parameter_fs(cp..., level)

		end


	return ( combine_functions_addrem(new_rmvs, rmv_internal_key),
					 combine_functions_addrem(new_adds, add_internal_param)
					 )
	

end 












#############################################################################
#module FileNames
#############################################################################

#usedkeys = ...
#DATAROOT defined in Constants 
#
#function get_FNG(input_dict::AbstractDict)::FilenameGenerator

#FilenameGenerator(usedkeys, input_dict, DATAROOT)
#end 

struct FilenameGenerator 

	# --- storage, for future function generation --- # 
	
	usedkeys::Union{<:Nothing, <:Vector{<:Symbol}, <:Function}

	digits::Union{<:Nothing, <:ODict}


	# --- actually used functions --- #  
	
	path::String  

	params_digits::Function 

	get_fname::Function

end 


get1(M::Module, k::Symbol) = Utils.getprop(M, k, nothing)

get1(M, k::Symbol) = M 

get1(M, ks...) = [get1(M, k) for k in ks]


# --- All 4 args are given --- # 

function FilenameGenerator(usedkeys::Union{<:Nothing,
																					 <:Vector{<:Symbol},
																					 <:Function}, 
													 digits::Union{<:Nothing, <:ODict},
													 path::AbstractString,
													 params_digits::Function,
													 )::FilenameGenerator

	function get_fname(args::UODicts)::Function
		
		fname(prefix([path, tostr(params_digits, args...)]))

	end  

	return FilenameGenerator(usedkeys, digits, prefix(path), 
													 params_digits, get_fname)

end 


# --- Only 3 args are given --- # 

function FilenameGenerator_(usedkeys::Union{<:Vector{<:Symbol}, <:Function}, 
													 digits::ODict,
													 path::Union{<:AbstractString,<:Tuple},
													 )::FilenameGenerator 

	FilenameGenerator(usedkeys, digits, prefix(path), 
										typical_params_digits(usedkeys, digits))

end 


function FilenameGenerator_(usedkeys::Union{<:Vector{<:Symbol}, <:Function}, 
													 digits::ODict,
													 params_digits::Function,
													 )::FilenameGenerator 

	FilenameGenerator(usedkeys, digits, "", params_digits)

end 


function FilenameGenerator_(usedkeys::Union{<:Vector{<:Symbol}, <:Function}, 
													 path::Union{<:AbstractString,<:Tuple},
													 params_digits::Function,
													 )::FilenameGenerator 

	FilenameGenerator(usedkeys, nothing, prefix(path), params_digits)


end 

function FilenameGenerator_(digits::ODict,
													 path::Union{<:AbstractString,<:Tuple},
													 params_digits::Function,
													 )::FilenameGenerator 

	FilenameGenerator(nothing, digits, prefix(path), params_digits)

end 




function FilenameGenerator(args::Vararg{<:Any,3})::FilenameGenerator

	function promising_candidate((FG,args), candidate::ODict)::Bool

		allunique(keys(candidate)) && length(candidate)<=length(args)
		
	end 
	
	
	function accept_sol((FG,args,), candidate::ODict)::Bool
	
		length(candidate)==length(args) || return false 
		
		return hasmethod(FG, typeof.(values(candidate)))
	
#usedkeys::Union{<:Nothing,
#																					 <:Vector{<:Symbol},
#																					 <:Function}, 
#													 digits::Union{<:Nothing, <:ODict},
#													 path::AbstractString,
#													 params_digits::Function,
	end 
	
	
	function possible_extensions((FG, args,), candidate::ODict)
	
		i = length(candidate)+1  

		i>length(args) && return []
	
		A = args[i]


		isa(A, Module) || return [merge(OrderedDict(), candidate, OrderedDict(i=>A)),]
	
		K = [(:usedkeys, :digits), (:digits, :path), (:path, :params_digits)][i]

		return map(intersect(names(A, all=true), setdiff(K, keys(candidate)))) do k
	
			merge(OrderedDict(), candidate, OrderedDict(k=>getproperty(A, k)))
	
		end 
	
	end 


	function output((FG, args), solutions, candidate)
	
#		push!(solutions, FG(values(candidate)...)) 
		push!(solutions, candidate)

		length(solutions)>1 && @warn "Ambiguous args."


#		if length(solutions)>1 
#			
#			println() 
#
#			for s in solutions
#
#				for (k,v) in s
#
#					println(k,"\t",v)
#
#				end
#				println()
#	
#			end 
#
#			error("\nAmbiguous args")
#		end 
	
	end 
	
	sols = Utils.Backtracking((FilenameGenerator_, args),
			 										  possible_extensions,
			 										  promising_candidate,
			 										  accept_sol,
			 										  output,
			 										  data->OrderedDict()
			 										  ) 

	length(sols)>1 && @warn "Random sol used"

#	return sols[1]
	
	return FilenameGenerator_(values(rand(sols))...)
	
end 




#
## --- Minimal input --- # 

function FilenameGenerator(path::Union{<:AbstractString, <:Module}, 
													 params_digits::Union{<:Function, <:Module}
													)::FilenameGenerator

	FilenameGenerator(nothing, nothing, 
										get1(path, :path),
										get1(params_digits, :params_digits))

end 




									
#import ...Utils

#jconst ROOT = "Data"


#using OrderedCollections: OrderedDict 

#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

# fname("sometask",(libA,pA1,pA2),(libB,pB),(libC,pC1,pC2,pC3,pC4))
# path = "root/sometask/signature(pA1)/signature(pA2)/sinature(pB)/..."


#function fname(name::AbstractString, func::Function, arg; kwargs...)::Function
#	strings = Utils.flatmapif(params_to_string, !isempty, [func(args)])
#
#	fname 
#	prefix = mkpath(join([ROOT;name; strings],"/"))
#
#	return (x::AbstractString="") -> string(prefix, "/", x) 
#
#end
#
#




tostr(x::AbstractString)::String = x

tostr(x::Symbol)::String = string(x)

tostr(x::Utils.List)::String = join(Utils.mapif(tostr, !isempty, x),"/")

tostr(M::Module)::String = tostr(filter(!isequal(:Main), Base.fullname(M)))



#function prefix(path::Union{<:AbstractString,<:Utils.List}, func::Function, 
#							 kwargs...)::String where N 
#
#	prefix(path; kwargs...)
#
#end 


function tostr(func::Function, arg::UODict)::String

	params_to_string(func(arg))

end 


function tostr(func::Function, args::UODicts; kwargs...)::String

	args  .|> func .|> params_to_string |> tostr

end 



#function fname(name::Union{AbstractVector,Tuple}; kwargs...)::Function
#
#	(args...; kw1...)-> fname(join(name,"/"), args...; kwargs..., kw1...)
#
#end

#function fname(M::Module, n::Int64=1, args...; kwargs...)::Function
#
#	fname(M, n, M.params_digits, args...; kwargs...)
#
#end
#
#
#function fname(M::Module, n::Int64, params_digits::Function, args...; kwargs...)::Function
#
#	fname(Base.fullname(M)[end-n+1:end], params_digits, args...; kwargs...)
#
#end


prefix(path)::String  = mkpath(tostr(path))
#function prefix(path)::String 
#	
#	s = tostr(path)
#
#	@show s 
#
#
#	@show pwd() 
#
#	@show ispath(pwd()*"/"*s)
#
#	mkpath(s) 
#
#	@show ispath(pwd()*"/"*s)
#	println()
#
#	return s 
#
#end 

prefix(path...)::String = prefix(path)



function fname(path::AbstractString)::Function 

	function strfname(x::AbstractString="")::String 
		
		tostr([path, x])

	end 

end 



#function fname(args...)::Function

#	prefix(args...)
							 

#
#===========================================================================#
#
# Turn one set of parameters into a string
#
#---------------------------------------------------------------------------#


function params_to_string(args::Vararg;
													separator::Union{<:AbstractString,<:Char}='-'
													)::String
	
	params, usedkeys, digits = params_digits_(args...)

	K = intersect(keys(digits), intersect(keys(params), usedkeys))

	isempty(K) && return ""

  return join([Utils.nr2string(params[k], digits[k]) for k in K], separator)

end

function params_to_string(args::Tuple; kwargs...)::String 

	params_to_string(args...; kwargs...)

end 








#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


get_usedkeys(usedkeys::Nothing, args...)::Nothing = nothing 

function get_usedkeys(usedkeys::AbstractVector{<:Symbol}, 
											args...)::Vector{Symbol}

	usedkeys

end 


get_usedkeys(usedkeys::Function)::Function = usedkeys

get_usedkeys(usedkeys::Function, P::UODicts)::Vector{Symbol} = usedkeys(P...)

function get_usedkeys(M::Module, args...)

	get_usedkeys(Utils.getprop(M, :usedkeys, Symbol[]), args...)

end 




#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#




function params_digits_(P::Tp, 
												uk::Union{<:AbstractVector{<:Symbol}, <:Function}, 
												digits::Td
											 )::Tuple{Tp, Vector{Symbol}, Td} where {Tp<:UODict, 
																															 Td<:ODict}
	(P, get_usedkeys(uk, P), digits)

end 


#function #params_digits_(P::Tp,
#					#							uk::Union{<:AbstractVector{<:Symbol}, <:Function}, 
#					#							digits::Td
#					#						 )::Tuple{Tp, Vector{Symbol}, Td} where {
#					#											Tp<:Tuple{UODicts},
#					#											Td<:ODict}
#
#	(P, get_#usedkeys(uk, P...), digits)
#
#end 




function params_digits_many(P::Utils.List, args...)::Tuple{<:UODict,
																													 Vector{Symbol},
																													 <:ODict}
	Utils.isList(P, UODict) || error()

	return params_digits_many(merge(P...), args...)

end 




function params_digits_many(P::Tp, args...
														)::Tuple{Tp, Vector{Symbol}, <:ODict
																		 } where Tp<:UODict 

	P_, usedkeys, digits = Utils.zipmap(a->get_params_digits(a,P), args)

	return params_digits_(P, union_usedkeys(usedkeys...), merge(digits...))

end






function typical_params_digits(usedkeys::Union{<:AbstractVector{Symbol},
																							 <:Function}=Symbol[],
															 digits::ODict=OrderedDict())::Function

	(P::UODict) -> params_digits_(P, usedkeys, digits)

end



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function union_usedkeys(uks_...)

	union_usedkeys_(Utils.mapif(get_usedkeys, !isnothing, uks_)...)

end 


function union_usedkeys_(uks...)

	isempty(uks) && return Symbol[]

	Utils.isList(uks, AbstractVector{<:Symbol}) && return union(uks...)

	if Utils.isList(uks, Union{<:Function, AbstractVector{<:Symbol}})
	
		function usedkeys(P::UODicts)::Vector{Symbol}

			union_usedkeys_((get_usedkeys(uk, P...) for uk in uks)...)

		end


		allusedkeys = union_usedkeys_(map(uks) do uk 

						u = get_usedkeys(uk)

						u isa AbstractVector{<:Symbol} && return u 

						u isa Function && return u()

						error("Wrong input")


		end...)


		usedkeys()::Vector{Symbol} = allusedkeys 



		return usedkeys


	end 

	error() 

end 



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function common_NrParamsSets(Ms::Vararg{<:Module})::Union{<:Nothing,<:Int}

	common_NrParamsSets([m for m in Ms])

end 



function common_NrParamsSets(Ms::AbstractVector{<:Module}
														)::Union{<:Nothing,<:Int}

	Ns = filter(>(0), Utils.getprop.(Ms, :NrParamSets, -1))

	isempty(Ns) && return nothing 

	N = Ns[1] 

	all(isequal(N), Ns) && return N 

	error()

end 



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


get_params_digits(F::Function)::Function = F

get_params_digits(F::Function, P::UODict)::Tuple = F(P)

get_params_digits(T::Tuple, args...)::Tuple = T






#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#




function combine_params_digits(args...)
	
#	N==1 && return merge_params_digits(args...)
				# apparently one argument is given, must be in fact a group
				
#	(P, get_usedkeys(uk, P), digits) 

	if all(args) do a 

		isa(a,ParamFlow) || isa(a,FilenameGenerator) || return false 

		return !isnothing(a.digits) && !isnothing(a.usedkeys)
			
		end 

		return typical_params_digits(union_usedkeys(args...),
																 merge((a.digits for a in args)...)
																 )
	end 


	tups = filter(a->isa(a,Tuple), args)

	!isempty(tups) && return params_digits_many(first.(tups), args...)

	return (P::UODict)->params_digits_many(P, args...)

end 


#
#function combine_params_digits(args::Vararg{<:Union{<:ParamFlow,<:FilenameGenerator}})
#	
#
#
#
#	function result(P...) 
#
#		tuples_ = map(args) do arg 
#			
#			typeof(arg)<:Module && return arg.params_digits(P...)
#			
#			typeof(arg)<:Function && return arg(P...)
#			
#			return arg
#			
#		end
#
#
#
#		out(tuples) = Tuple(map(enumerate([merge, union, merge])) do (i,f)
#	
#													f([t[i] for t in tuples]...)
#
#												end)
#
#
##		tuples_ = ((p1,u1,d1), (p2,u2,d2) etc)
#
#		tuples_[1][1] isa UODict && return out(tuples_)
#
#	
##		tuples_ = ( [ (pi1,ui1,di1) etc ], [ (pj1,uj2,dj2) etc ], etc  )
#
#		return [out(tuples) for tuples in zip(tuples_...)]
#
##		out( (pi1,ui1,di1), (pj1,uj1,dj1) ), out( (pi2,ui2,...), (j2) )
#
#	end
#
#
#
#	for arg in args
#
#		typeof(arg)<:Module && continue
#
#		typeof(arg)<:Function && continue
#
#		(p, u, d) = arg
#		
#		return result(p)
#
##		if at least one tuple (p, u ,d) is given, will return the total tuple 
##								(total_p, total_u, total_d)
#			
#	end
#
#	return (P...) -> result(P...)
#	
#
#	
## if only modules are given, 
##						 return the function which computes the total tuple
#
#end
#
##
##
##
##
#
#
#
#
#
#############################################################################
#end
#############################################################################

struct ParamFlow

	# --- own properties --- # 

	NrParamSets::Union{<:Nothing, <:Int}

	allparams::Function



	# --- storage, for future function generation --- # 
	
	usedkeys::Union{<:Nothing, <:Vector{<:Symbol}, <:Function}

	digits::Union{<:Nothing, <:ODict}


	# --- actually used functions --- #  
	
	path::String  

	params_digits::Function 

	get_fname::Function 



	# --- #



	function ParamFlow(N::Union{<:Nothing,<:Int}, allparams::Function, 
										 FNG::FilenameGenerator)::ParamFlow

		new(N, allparams, (getproperty(FNG, p) for p in propertynames(FNG))...)
	
	end 



	function ParamFlow(allparams::Function, FNG::FilenameGenerator)::ParamFlow

		ParamFlow(nothing, allparams, FNG)

	end 



	function ParamFlow(allparams::Function, arg1, arg2, args...)::ParamFlow
	
		ParamFlow(allparams, FilenameGenerator(arg1, arg2, args...))
	 
	end


	function ParamFlow(N::Union{<:Nothing, <:Int}, allparams::Function, 
										 arg1, arg2, args...)::ParamFlow

		ParamFlow(N, allparams, FilenameGenerator(arg1, arg2, args...))

	end 




	function ParamFlow(NrParamSets::Int, input_dict::UODict, args...)::ParamFlow
	
		FNG = FilenameGenerator(args...) 
	
		function allparams(args::Vararg{<:UODict,N}) where N 
	
			N < NrParamSets || error()  

			return typical_allparams(FNG.usedkeys, input_dict)
	
		end 

		return ParamFlow(NrParamSets, allparams, FNG)
	 
	end



	function ParamFlow(M::Module, input_dict::UODict, args...)::ParamFlow

		ParamFlow(M.NrParamSets, input_dict, args...)

	end 


	function ParamFlow(input_dict::UODict, args...)::ParamFlow 

		ParamFlow(1, input_dict, args...)

	end 






#	function ParamFlow(N::Union{<:Nothing,<:Int},
#										 allparams::Function, arg1, arg2, args...)::ParamFlow
#
#		ParamFlow(N, allparams, FilenameGenerator(arg1, arg2, args...))
#
#	end 

end 

#usedkeys = ...
#NrParamSets = ...
#
#DATAROOT defined in Constants 
#
#function get_PF(input_dict::AbstractDict)::FilenameGenerator

#ParamFlow((NrParamSets, input_dict), usedkeys, input_dict, DATAROOT)
#end 



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#





function get_usedkeys(S::Union{<:FilenameGenerator, <:ParamFlow}, args...)

	get_usedkeys(S.usedkeys, args...)

end 


function get_params_digits(S::Union{<:FilenameGenerator, <:ParamFlow}, 
													 args...)

	get_params_digits(S.params_digits, args...)

end 




#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function good_methods(PF::ParamFlow, F::Function, args...)::Bool
	
	for i in 1:Utils.Assign_Value(PF.NrParamSets, 1)

		hasmethod(F, NTuple{i,<:AbstractDict}, args...) || return false 

	end 

	return true 

end 



function add_kwargs(PF::ParamFlow, F::Function; kwargs...)::Function

	new_kwargs = [:get_fname, ]

	@assert good_methods(PF, F) 


	isempty(new_kwargs) && isempty(kwargs) && return F 


	ks1 = filter(propertynames(PF)) do k 

		in(k, new_kwargs) && good_methods(PF, F, tuple(k)) 

	end 


	ks2 = filter(collect(keys(kwargs))) do k

		good_methods(PF, F, tuple(k))

	end 

	isempty(ks1) && isempty(ks2) && return F  # it accepts no kwargs


	kwa1 = NamedTuple{ks1}(getproperty.([PF], ks1)) 

	kwa2 = NamedTuple{Tuple(ks2)}(kwargs[k] for k in ks2)


	return function f(args::Vararg{<:UODict}; kwa...)
	
		F(args...; kwa1..., kwa2..., kwa...)

	end 

end 




struct Calculation

	name::String

	PF::Union{<:ParamFlow,<:FilenameGenerator}
	
	Compute::Function 

	Read::Function 

	FoundFiles::Function 


	function Calculation(name::Union{<:AbstractString,<:Symbol},
											 PF::Union{<:ParamFlow, <:FilenameGenerator},
											 CRF::Vararg{<:Function,3}; kwargs...)

		new(string(name), 
				PF, (add_kwargs(PF, f; kwargs...) for f in CRF)...)

	end 


	function Calculation(name::Union{<:AbstractString,<:Symbol},
											 PF::Union{<:ParamFlow, <:FilenameGenerator},
											 Compute::Function; kwargs...)

		C = add_kwargs(PF, Compute; kwargs...) 

		return new(string(name), PF, C, C, (args...)->true)

	end 



	function Calculation(name::Union{<:AbstractString,<:Symbol},
											 PF::Union{<:ParamFlow, <:FilenameGenerator},
											 Compute::Function, ::Nothing, ::Any; kwargs...)

		Calculation(name, PF, Compute; kwargs...) 

	end 


	function Calculation(name::Union{<:AbstractString,<:Symbol},
											 PF::Union{<:ParamFlow, <:FilenameGenerator},
											 Compute::Function, ::Nothing, ::Nothing; kwargs...)

		Calculation(name, PF, Compute; kwargs...)

	end 


	function Calculation(PF::Union{<:ParamFlow, <:FilenameGenerator},
											 Compute::Function, args...; kwargs...)

		Calculation("", PF, Compute, args...; kwargs...)

	end 




	function Calculation(name::Union{<:AbstractString,<:Symbol},
											 PF::Union{<:ParamFlow, <:FilenameGenerator},
											 M::Module; kwargs...)

		Calculation(name, PF, Utils.getprop(M,[:Compute,:Read,:FoundFiles])...;
								kwargs...)

	end 


	function Calculation(PF::Union{<:ParamFlow, <:FilenameGenerator},
											 M::Module; kwargs...)

		Calculation(tostr(M), PF, M; kwargs...)

	end 





#	function Calculation(name::AbstractString, args...)


#	end 


end 



tostr(C::Calculation)::String = C.name 

#===========================================================================#
#
#	Depending on the kwargs, converts the parameter set into a flat list, 
#			 the other way around, or fills a missing parameter
#	
#		 kwargs must contain:
#		 	- M = the module for which the calculation is done
#		 	- get_new = gives the parameters at a certain level
#		 	- convert_one = method to convert one dictionary
#		 	- convert_many = specifies how to distribute 'convert_one' on the
#		 										items in a dictionary list
#		 	- repl = can replace some parameter with a given value (to fix it)
#		 	- final_f = to obtain the result in a specific manner
#
#---------------------------------------------------------------------------#

#			M Must have M.NrParamSets M.allparams

get_NrPSets(m::Module)::Int = m.NrParamSets
get_NrPSets(pf::ParamFlow)::Int = pf.NrParamSets 
get_NrPSets(c::Calculation)::Int = get_NrPSets(c.PF)


get_allP(m::Module, args...) = m.allparams(args...) 
get_allP(pf::ParamFlow, args...) = pf.allparams(args...) 
get_allP(m::Calculation, args...) = get_allP(m.PF, args...)



function convert_params(M::Union{<:Module, <:Calculation},
												args=[];
												get_new=(l,a)->get_allP(M, a...),
												repl=nothing,
												convert_one=identity,
												convert_many=map,
												final_f=identity)

	function act_on_plev(f1::T1) where T1<:Function
	
		return function (arg::T2) where T2
	
			T2 <: AbstractDict && return f1(arg)
		
			Utils.isList(T2,AbstractDict) && return convert_many(f1,arg)
		
			error("\n $arg \nParameter structure not understood: $T2")
	
		end
	
	end


	

	for level in 1:get_NrPSets(M)

		new = get_new(level,args)

		if isnothing(repl) push!(args,new) else 

				push!(args, act_on_plev(a->repl(level,a))(new))

		end

	end


	return final_f(map(act_on_plev(convert_one),args))

end

#===========================================================================#
#
# Particularizes convert_params in order to obtain a single dictionary
# 		with unique, string keys, from a structured set P
#
# 		(the inverse operation of convertParams_fromPlot)
#
#---------------------------------------------------------------------------#

function convertParams_toPlot(M::Union{<:Module, <:Calculation},
															P::Utils.List; kwargs...)

	isempty(P) && return convertParams_toPlot(M; kwargs...)

	return convertParams_toPlot(M; get_new = (level,args)->P[level], kwargs...)

end 


function convertParams_toPlot(M::Union{<:Module, <:Calculation},
															P::Nothing=nothing; 
															kwargs...)::OrderedDict

	convert_params(M;

		convert_one = convertPairs_toPlot,

		convert_many = Utils.flatmap,

		final_f = args -> OrderedDict(vcat(args...)),

		kwargs...

		) 

end


#===========================================================================#
#
#	Particularizes convert_params in order to obtain a structured
#			parameter set from a single, flat, dictionary P
#
#			(the inverse operation of convertParams_toPlot)
#
#---------------------------------------------------------------------------#

function convertParams_fromPlot(M::Union{<:Module,<:Calculation},
																P::UODict; kwargs...)

	convert_params(M;

		repl = function (level::Int, params::UODict)

										K = Utils.mapif(k1->(k1,convertKey_toPlot(params,k1)),

																		k -> haskey(P,k[2]),

																		collect(keys(params)))

										return OrderedDict(k1 => P[k2] for (k1,k2) in K)

									end,

		kwargs...

		)

end


##################### REFORMULTE WITH Utils.Backtracking !! 

#===========================================================================#
#
#	Fill internal parameters
#
#---------------------------------------------------------------------------#


function fillParams_internalP(M::Union{<:Module, <:Calculation},
															P::Utils.List, ::Nothing;
															kwargs...)
	P 

end 


function fillParams_internalP(M::Union{<:Module, <:Calculation},
															P::Utils.List, add_internal_param::Function;
															kwargs...)

	convert_params(M;
								 
								 get_new = (level,args) -> P[level],

								 repl = add_internal_param,

								 kwargs...
	
								)
end 



function fillParams_internalPs(M::Union{<:Module, <:Calculation}, 
															 P::Utils.List, ::Nothing;
															kwargs...)
	P

end 


function fillParams_internalPs(M::Union{<:Module, <:Calculation}, 
															 P::Utils.List, 
															 add_internal_params::Utils.List;
															 kwargs...)

	map(add_internal_params) do add

					fillParams_internalP(M, P, add; kwargs...)

	end

end


#function f_fillParams_internal(M::Module, add_internal_param, add_internal_params)
#
#	function (P; shuffle=false)
#
#		P1 = fillParams_internalP(M, P, add_internal_param)
#
#		shuffle = !isnothing(add_internal_params) & shuffle
#
#		return fillParams_internalPs(M, P1, 
#								(shuffle ? Random.shuffle : identity)(add_internal_params)
#															)
#	end 
#end 


#===========================================================================#
#
#	From module M, get parameter combinations: 
#									- all (if kwargs not given)
#									- with contraint/condition 'cond'
#									- with some fixed parameter, implemented by 'repl' 
#
#---------------------------------------------------------------------------#

function get_paramcombs(M::Union{<:Module, <:ParamFlow, <:Calculation};
												cond=nothing, repl=nothing) # First call

	get_paramcombs(M, 1, [];
								 cond = combine_functions_cond(cond),
								 repl = combine_functions_addrem(repl),
								)
end 



function get_paramcombs(M::Union{<:Module, <:ParamFlow, <:Calculation},
												level::Int, old::AbstractVector; 
												cond::Function, repl::Function) # Further calls 

	#M must have M.allparams::Function, M.NrParamSets::Int,

	news = Utils.mapif(pcomb -> vcat(old..., pcomb),

										 new->cond(new, level),

										 Utils.AllValCombs(repl(level, get_allP(M, old...))))


	level == get_NrPSets(M) && return news

	return Utils.flatmap(news) do new 
				
						get_paramcombs(M, level+1, new; cond=cond, repl=repl) 
	end

end




function f_get_plotparams(M::Union{<:Module, <:Calculation}, 
													rmv_internal_key::Nothing=nothing)

#	getpp() = convertParams_toPlot(M) 

	getpp(P::UODicts) = convertParams_toPlot(M, P)
												
#	return getpp 

end 


function f_get_plotparams(M::Union{<:Module, <:Calculation}, 
													rmv_internal_key::Function)

	getpp(P::UODicts) = convertParams_toPlot(M, P; 
																					 repl=rmv_internal_key) 

end 



function f_fill_internal(M::Union{<:Module, <:Calculation}, args...)
	
	fill_int(P::Utils.List) = Parameters.fillParams_internalP(M, P, args...)
	
end 


#############################################################################
end
