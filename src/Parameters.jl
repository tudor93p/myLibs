module Parameters
#############################################################################

import ...OrderedDict, ...Utils 
import Random

export UODict, ODict

const ODict = Union{<:OrderedDict, <:NamedTuple}
const UDict = AbstractDict
const UODict = Union{<:AbstractDict, <:OrderedDict, <:NamedTuple}
#const UODicts = Vararg{<:UODict}


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





#===========================================================================#
#
# applies several add or remove functions recursively
#
#---------------------------------------------------------------------------#


function apply_one_or_many(f::Function, l::Int, P...; default)

	applicable(f, l, P...) && return f(l, P...)

	applicable(f, l, P[l]) && return f(l, P[l])

	@warn "Function $f not applicable. Returning P[l]"

	return default 

	error()

end 


function recrepl(l::Int, P...; fs::AbstractVector{Function}=Function[]) 

	isempty(fs) && return P[l]

	new_Pl = apply_one_or_many(fs[1], l, P...; default=P[l])

	return recrepl(l, (i==l ? new_Pl : p for (i,p) in enumerate(P))...;
								 fs=fs[2:end]) 

end 



function recrepl(fs::AbstractVector)::Function

	recrepl_(l::Int, P...) = recrepl(l, P...; fs=Vector{Function}(fs))

end 



function combine_functions_addrem(args...)::Function

	recrepl(Utils.flat(args...; keep=x->isa(x,Function)))

end 




#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#



function combine_functions_cond(args...)::Function

	good_fs = Utils.flat(args...; keep=x->isa(x,Function))

	function combined_cond(l::Int, P...)::Bool

		for cond in good_fs
			
			apply_one_or_many(cond, l, P...; default=true) || return false 

		end 

		return true 

	end


#	function combined_cond(P::Utils.List, l::Int)::Bool 
#
#		all(combined_cond.(P,l))
#
#	end 


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

function convertKey_toPlot(k; kwargs...)::String 

	string(k)

end

function convertKey_toPlot(p::UODict, k; separator="__")::String 

	haskey(p,:Label) && return join([p[:Label], k], separator) 

	return convertKey_toPlot(k; separator=separator)

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


function replace_parameter_fs(getP, Keys, levels)::Tuple{Function,Function}

	replace_parameter_fs_(getP, Utils.flat(Keys), in(Utils.flat(levels)))

end



function replace_parameter_fs_(getP::Function, Keys::AbstractVector, 
															 act::Function 
															 )::Tuple{Function,Function}

	rmv,aux = replace_parameter_fs_(nothing, Keys, act) 

	m(p::UODict, q::Utils.List) = Utils.adapt_merge(p, zip(Keys,q))

	m(p::Utils.List, q) = m.(p, q)

	add(l::Int, P...) = act(l) ? m(P[l], getP(l, P...)) : P[l]
	
	return rmv,add

end 



function replace_parameter_fs_(getP::Nothing,
															 Keys::AbstractVector, 
															 act::Function
															)::Tuple{Function,Function}

	f(p::UODict) = Utils.dict_diff(p, Keys)

	f(P::Utils.List) = f.(P)


	rmv(l::Int, P...) = act(l) ? f(P[l]) : P[l]

	add(l::Int, P...) = P[l]

	return rmv,add

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

	if isnothing(constrained_params) || Utils.isListLen(constrained_params, 0)

		return ( combine_functions_addrem(rmv_internal_key),
						 combine_functions_addrem(add_internal_param)
						 ) 
	end 


	new_rmvs,new_adds = Utils.zipmap(constrained_params) do (level,cp)
			
				isempty(cp) ? (nothing,nothing) : replace_parameter_fs(cp..., level)

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

#struct FilenameGenerator # will be deleted

#
#	# --- storage, for future function generation --- # 
#	
#	usedkeys::Union{<:Nothing, <:Vector{<:Symbol}, <:Function}
#
#	digits::Union{<:Nothing, <:ODict}
#
#
#	# --- actually used functions --- #  
#	
#	path::String  
#
#	params_digits::Function 
#
#	get_fname::Function
#
#end 



#
#function FilenameGenerator(usedkeys::Union{<:Nothing, <:Vector{<:Symbol}, <:Function},
#													 digits::Union{<:Nothing, <:ODict},
#													 path::String,
#													 params_digits::Function,
#													 get_fname::Function
#													 )
#
#	usedkeys, digits, path, params_digits, get_fname
#
#end 
#
#
#get1(M::Module, k::Symbol) = Utils.getprop(M, k, nothing)
#
#get1(M, k::Symbol) = M 
#
#get1(M, ks...) = [get1(M, k) for k in ks]
#
#
## --- All 4 args are given --- # 
#
#function FilenameGenerator(usedkeys::Union{<:Nothing,
#																					 <:Vector{<:Symbol},
#																					 <:Function}, 
#													 digits::Union{<:Nothing, <:ODict},
#													 path::AbstractString,
#													 params_digits::Function,
#													 )#::FilenameGenerator



function construct_get_fname(path_, 
														 all_params_digits::Vararg{<:Function}
														 )::Function 

	path = prefix(path_) 

	return function get_fname(args...)::Function

		pref = prefix(path, zip(all_params_digits, Utils.flat(args))...)

#		pref = prefix(path, zip(all_params_digits, args)...)

		return function strfname(x::AbstractString="")::String 
			
			tostr([pref, x])
	
		end 
	
	end 

end 


function construct_get_fname(path,#::AbstractString,
														 ingreds_pds...
														 )::Function 

	construct_get_fname(path,
#											vcat(typical_params_digits(ingreds_pds...))...)
											Utils.flat(typical_params_digits(ingreds_pds...))...)


end 






#	fname(prefix([path, tostr(params_digits, args...)])) 




#	return FilenameGenerator(usedkeys, digits, prefix(path), 
#													 params_digits, get_fname)
#
#end 
#
#
## --- Only 3 args are given --- # 
#
#function FilenameGenerator_(usedkeys::Union{<:Vector{<:Symbol}, <:Function}, 
#													 digits::ODict,
#													 path::Union{<:AbstractString,<:Tuple},
#													 )#::FilenameGenerator 
#
#	FilenameGenerator(usedkeys, digits, prefix(path), 
#										typical_params_digits(usedkeys, digits))
#
#end 
#
#
#function FilenameGenerator_(usedkeys::Union{<:Vector{<:Symbol}, <:Function}, 
#													 digits::ODict,
#													 params_digits::Function,
#													 )#::FilenameGenerator 
#
#	FilenameGenerator(usedkeys, digits, "", params_digits)
#
#end 
#
#
#function FilenameGenerator_(usedkeys::Union{<:Vector{<:Symbol}, <:Function}, 
#													 path::Union{<:AbstractString,<:Tuple},
#													 params_digits::Function,
#													 )#::FilenameGenerator 
#
#	FilenameGenerator(usedkeys, nothing, prefix(path), params_digits)
#
#
#end 
#
#function FilenameGenerator_(digits::ODict,
#													 path::Union{<:AbstractString,<:Tuple},
#													 params_digits::Function,
#													 )#::FilenameGenerator 
#
#	FilenameGenerator(nothing, digits, prefix(path), params_digits)
#
#end 
#
#
#
#
#
#
#
#function FilenameGenerator(args::Vararg{<:Any,3})#::FilenameGenerator
#
#	function promising_candidate((FG,args), candidate::ODict)::Bool
#
#		allunique(keys(candidate)) && length(candidate)<=length(args)
#		
#	end 
#	
#	
#	function accept_sol((FG,args,), candidate::ODict)::Bool
#	
#		length(candidate)==length(args) || return false 
#		
#		return hasmethod(FG, typeof.(values(candidate)))
#	
#	end 
#	
#	
#	function possible_extensions((FG, args,), candidate::ODict)
#	
#		i = length(candidate)+1  
#
#		i>length(args) && return []
#	
#		A = args[i]
#
#
#		isa(A, Module) || return [merge(OrderedDict(), candidate, OrderedDict(i=>A)),]
#	
#		K = [(:usedkeys, :digits), (:digits, :path), (:path, :params_digits)][i]
#
#		return map(intersect(names(A, all=true), setdiff(K, keys(candidate)))) do k
#	
#			merge(OrderedDict(), candidate, OrderedDict(k=>getproperty(A, k)))
#	
#		end 
#	
#	end 
#
#
#	function output((FG, args), solutions, candidate)
#	
##		push!(solutions, FG(values(candidate)...)) 
#		push!(solutions, candidate)
#
#		length(solutions)>1 && @warn "Ambiguous args."
#
#	
#	end 
#	
#	sols = Utils.Backtracking((FilenameGenerator_, args),
#			 										  possible_extensions,
#			 										  promising_candidate,
#			 										  accept_sol,
#			 										  output,
#			 										  data->OrderedDict()
#			 										  ) 
#
#	length(sols)>1 && @warn "Random sol used"
#
##	return sols[1]
#	
#	return FilenameGenerator_(values(rand(sols))...)
#	
#end 
#
#
#
#
##
### --- Minimal input --- # 
#
#function FilenameGenerator(path::Union{<:AbstractString, <:Module}, 
#													 params_digits::Union{<:Function, <:Module}
#													)#::FilenameGenerator
#
#	FilenameGenerator(nothing, nothing, 
#										get1(path, :path),
#										get1(params_digits, :params_digits))
#
#end 
#
#
#
#
									
#import ...Utils

#jconst ROOT = "Data"


#using OrderedCollections: OrderedDict 


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

function tostr(x::Utils.List)::String  

	s = Utils.mapif(!isempty, collect(x)) do xi 

			Utils.isList(xi) ? tostr(xi...) : tostr(xi) 

	end 

	return join(s,"/")

end 

tostr(M::Module)::String = tostr(filter(!isequal(:Main), Base.fullname(M)))


function tostr(func::Function, arg::UODict)::String

	params_to_string(func(arg))

end 


#function tostr(func::Function, args::UODicts; kwargs...)::String
#
#	args  .|> func .|> params_to_string |> tostr
#
#end 



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

prefix(path...)::String = prefix(path)




#===========================================================================#
#
# Turn one set of parameters into a string
#
#---------------------------------------------------------------------------#


function params_to_string(args...;
													separator::Union{<:AbstractString,<:Char}='-'
													)::String
	
	params, usedkeys, digits = params_digits_(args...)

	K = collect(intersect(keys(digits), intersect(keys(params), usedkeys)))

	isempty(K) && return ""

#	return join(vcat((Utils.nr2string(params[k], digits[k]) for k in K)...), separator)

	return join(Utils.flatmap(k->Utils.nr2string(params[k], digits[k]), K), 
							separator)

end 



function params_to_string(args::Tuple; kwargs...)::String 

	params_to_string(args...; kwargs...)

end 








#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


#get_usedkeys(usedkeys::Nothing, args...)::Nothing = nothing 
#
function get_usedkeys(usedkeys::AbstractVector{<:Symbol}, 
											args...)::Vector{Symbol}

	usedkeys

end 


get_usedkeys(usedkeys::Function)::Function = usedkeys

get_usedkeys(usedkeys::Function, P::UODict)::Vector{Symbol} = usedkeys(P)

function get_usedkeys(M::Module, args...)
#
	get_usedkeys(Utils.getprop(M, :usedkeys, Symbol[]), args...)
#	get_usedkeys(M.usedkeys, args...)


end 




#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#




function params_digits_(P::Tp, 
												uk::Union{<:AbstractVector{<:Symbol}, 
																	<:Function,
#																	<:Module 
																	}, 
												digits::Td
											 )::Tuple{Tp, Vector{Symbol}, Td} where {Tp<:UODict, 
																															 Td<:ODict}
	(P, get_usedkeys(uk, P), digits)

end 




#function params_digits_many(P::Utils.List, args...)::Tuple{<:UODict,
#																													 Vector{Symbol},
#																													 <:ODict}
#	Utils.isList(P, UODict) || error()
#
#	return params_digits_many(merge(P...), args...)
#
#end 
#
#
#
#
#function params_digits_many(P::Tp, args...
#														)::Tuple{Tp, Vector{Symbol}, <:ODict
#																		 } where Tp<:UODict 
#
#	P_, usedkeys, digits = Utils.zipmap(a->get_params_digits(a,P), args)
#
#	return params_digits_(P, union_usedkeys(usedkeys...), merge(digits...))
#
#end
#
#
#

#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


#function typical_params_digits(M::Module, digits::ODict)::Function 




function typical_params_digits(usedkeys::Union{<:AbstractVector{Symbol},
																							 <:Function},
															 digits::ODict)::Function

	(P::UODict) -> params_digits_(P, usedkeys, digits)

end


function typical_params_digits(usedkeys::Union{<:AbstractVector{Symbol},
																							 <:Function},
															 digits::AbstractVector{<:ODict}
															 )::Vector{Function}

	[typical_params_digits(usedkeys, d) for d in digits]

end



function typical_params_digits(M::Module, digits::ODict)::Function

	typical_params_digits(M.usedkeys, digits)

end 



function typical_params_digits(M::Module)::Function 

	M.params_digits

end 



function typical_params_digits(pd::Function)::Function 

	pd 

end 



#function typical_params_digits(::Nothing, ::Any, pd::Function)::Function 
#
#	pd
#
#end



#function typical_params_digits(t::Tuple)::Function 
#
#	typical_params_digits(t...)
#
#end 
#

function typical_params_digits(tuples::Vararg{<:Tuple})::Vector{Function}

	Utils.flat([typical_params_digits(t...) for t in tuples])

end 


#M.usedkeys 
#M.params_digits 
#M.allparams 



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


#
#function construct_get_allparams(fs::Vararg{<:Function, NrArgs}
#																)::Function where NrArgs 
#
#
#end 
#
#
#function construct_get_allparams(NrParamSets::Int, 
#																 fs::Vararg{<:Function,NrArgs}
#																 )::Function where NrArgs
#	
#	
#	@assert NrParamSets==NrArgs
#
#	return construct_get_allparams(fs...)
#
#end 
#
#
#
#function construct_get_allparams(args...)::Function
#
#	construct_get_allparams(vcat(typical_allparams(args...))...)
#
#end 
#

#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function typical_allparams(allparams::Function)::Function
	
	allparams

end 

function typical_allparams(M::Module)::Function
	
	M.allparams

end 







function typical_allparams(NrParamSets::Int, allparams::Function)::Function

	@assert hasmethod(allparams, NTuple{NrParamSets-1, <:UODict})

	return typical_allparams(allparams)

end 



function keep_usedkeys(usedkeys::Function, allp::UODict)::UODict 

	keep_usedkeys(usedkeys(), allp)

end 

function keep_usedkeys(usedkeys::AbstractVector{<:Symbol},
																allp::UODict)::UODict 

	Utils.dict_keepkeys(allp, usedkeys)

end  



function typical_allparams(NrParamSets::Int,
													 usedkeys::AbstractVector{<:Symbol},
													 allp::UODict)::Function 

	d = keep_usedkeys(usedkeys, allp)

	return function allparams(P...)::UODict#Dict{Symbol,Any}# where {T,N}
	
		N = length(P)

		N==NrParamSets-1 && return d 
		
		error("The function can only be called with ",NrParamSets-1," arguments. You provided $N.")
		
	end 

end


function typical_allparams(NrParamSets::Int,
													 usedkeys::AbstractVector{<:Symbol},
													 AllP::AbstractVector{<:UODict})::Function 

	D = [keep_usedkeys(usedkeys, allp) for allp in AllP]

	return function allparams(P...)::Vector{UODict}

		N = length(P) 

		N==NrParamSets-1 && return D 
		
		error("The function can only be called with ",NrParamSets-1," arguments. You provided $N.")
		
	end 

end


function typical_allparams(NrParamSets::Int, 
													 usedkeys::Function, 
													 allp::Union{<:UODict, <:AbstractVector{<:UODict}}
													 )::Function

	typical_allparams(NrParamSets, usedkeys(), allp)

end 

function typical_allparams(NrParamSets::Int, 
													 M::Module, 
													 allp::Union{<:UODict, <:AbstractVector{<:UODict}}
													 )::Function

	typical_allparams(NrParamSets, M.usedkeys, allp)

end 


function typical_allparams(usedkeys::Union{<:AbstractVector{<:Symbol},
																					 <:Module,
																					 <:Function},
													 allp::Union{<:UODict, <:AbstractVector{<:UODict}},
													 )::Function 

	typical_allparams(1, usedkeys, allp)

end 







#typical_allparams(t::Tuple)::Function = typical_allparams(t...)

function typical_allparams(tuples::Vararg{<:Tuple, NF})::Function where NF

	fs = Dict{Int,Function}( 

						map(enumerate(tuples)) do (i,t)

							t[1] isa Int && return t[1] => typical_allparams(t...)

							return i => typical_allparams(i, t...)

						end)


	return function allparams(P...)
	
		N = length(P) 

		N<NF && return fs[N+1](P...)

		error("The function can be called with at most ",NF-1," arguments. You provided $N")

	end 

end



function typical_allparams(NrParamSets::Int,
													 tuples::Vararg{<:Tuple,NF})::Function where NF

	@assert NrParamSets==NF "The number of tuples provided ($NF) should match the claimed 'NrParamSets' ($NrParamSets)"

	return typical_allparams(tuples...)

end




#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#



function sep_pdap(usedkeys::Union{<:AbstractVector{Symbol}, <:Function},
									digits::Union{<:ODict,<:AbstractVector{<:ODict}},
									allp::Union{<:UODict,<:AbstractVector{<:UODict}},
									)::Tuple{Tuple,Tuple}

	(usedkeys, digits), (usedkeys, allp)

end 


function sep_pdap(usedkeys::Union{<:AbstractVector{Symbol}, <:Function},
									digits::Union{<:ODict,<:AbstractVector{<:ODict}},
									arg::Union{<:Function, <:Module})::Tuple{Tuple,Tuple}

	(usedkeys, digits), (arg,)

end 

function sep_pdap(pd::Union{<:Function, <:Module},
									usedkeys::Union{<:AbstractVector{Symbol}, <:Function},
									allp::Union{<:UODict,<:AbstractVector{<:UODict}},
									)::Tuple{Tuple,Tuple}

	(pd,),(usedkeys,allp)

end  


function sep_pdap(usedkeys::Module,
									digits::Union{<:ODict,<:AbstractVector{<:ODict}},
									allp::Function)::Tuple{Tuple,Tuple}

	(usedkeys, digits), (allp,)

end 


function sep_pdap(M::Module, 
									digits::Union{<:ODict,<:AbstractVector{<:ODict}},
									allp::Union{<:UODict,<:AbstractVector{<:UODict}},
									)::Tuple{Tuple,Tuple}

	(M, digits), (M, allp)

end 






#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function sep_pdap(M::Module)::Tuple{Tuple,Tuple}

	(M,),(M,)

end 

#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#



function sep_pdap(pd::Union{<:Function,<:Module},
									allp::Union{<:Function,<:Module})::Tuple{Tuple,Tuple}

	(pd,),(allp,)
	
end 



function sep_pdap(M::Module, 
									digits::Union{<:ODict,<:AbstractVector{<:ODict}},
									)::Tuple{Tuple,Tuple}

	(M, digits), (M,)

end 




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


function union_usedkeys(uks_...)

	union_usedkeys_(Utils.mapif(get_usedkeys, !isnothing, uks_)...)

end 

union_usedkeys_()::Vector{Symbol} = Symbol[]


function union_usedkeys_(uks...)

	Utils.isList(uks, AbstractVector{<:Symbol}) && return union(uks...)

	if Utils.isList(uks, Union{<:Function, AbstractVector{<:Symbol}})
	
		function usedkeys(P::UODict)::Vector{Symbol}

			union_usedkeys_((get_usedkeys(uk, P) for uk in uks)...)

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


#get_params_digits(F::Function)::Function = F
#
#get_params_digits(F::Function, P::UODict)::Tuple = F(P)
#
#get_params_digits(T::Tuple, args...)::Tuple = T






#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#




#function combine_params_digits(args...)
#	
##	N==1 && return merge_params_digits(args...)
#				# apparently one argument is given, must be in fact a group
#				
##	(P, get_usedkeys(uk, P), digits) 
#
#	if all(args) do a 
#
#		isa(a, ParamFlow) || isa(a,FilenameGenerator) || return false 
#
#		return !isnothing(a.digits) && !isnothing(a.usedkeys)
#			
#		end 
#
#		return typical_params_digits(union_usedkeys(args...),
#																 merge((a.digits for a in args)...)
#																 )
#	end 
#
#
#	tups = filter(a->isa(a,Tuple), args)
#
#	!isempty(tups) && return params_digits_many(first.(tups), args...)
#
#	return (P::UODict)->params_digits_many(P, args...)
#
#end 





























struct ParamFlow

	# --- own properties --- # 

	NrParamSets::Int


#	# --- storage, for future function generation --- # 
#	path::String  
#	
#	usedkeys::Union{<:Nothing, <:Vector{<:Symbol}, <:Function}
#
#	digits::Union{<:Nothing, <:ODict}

#	params_digits::Vector{Function}

	# --- actually used functions --- #  
	

	get_fname::Function  

	allparams::Function


	constraint::Union{<:Nothing, <:Function} 

	adjust::Union{<:Nothing, <:Function} 

#	adjust::Function

	# --- #


	function ParamFlow(path, NrParamSets::Int, 
										 tuples::Vararg{<:Tuple,N}; 
										 constraint=nothing,
										 adjust=nothing,
										 kwargs...) where N

		@assert NrParamSets==N "The number of tuples provided ($N) should match the claimed 'NrParamSets' ($NrParamSets)"

		args_pd, args_ap = zip((sep_pdap(t...) for t in tuples)...) 

		return new(NrParamSets,
							 construct_get_fname(path, args_pd...),
							 typical_allparams(NrParamSets, args_ap...),
							 combine_functions_cond(constraint),
							 combine_functions_addrem(adjust),
							 )

	end 


	function ParamFlow(path, tuples::Vararg{<:Tuple,N}; kwargs...) where N

		ParamFlow(path, N, tuples...; kwargs...)

	end 


	function ParamFlow(path, NrParamSets::Int, args...; 
										 constraint=nothing, adjust=nothing, kwargs...)

		@assert NrParamSets==1 "The arguments provided contradict the claimed 'NrParamSets' ($NrParamSets)"

		args_pd,args_ap = sep_pdap(args...)

		return new(NrParamSets,
							 construct_get_fname(path, args_pd...),
							 typical_allparams(NrParamSets, args_ap...),
							 combine_functions_cond(constraint),
							 combine_functions_addrem(adjust),
							 )

	end 

	function ParamFlow(path, args...; kwargs...)

		ParamFlow(path, 1, args...; kwargs...)

	end 

#beginning: path (and NrParamSets)
#

#next: tuples or args 
#

	
#	function ParamFlow(NrParamSets::Int, 
#										 path, 
#										 uk::Union{<:Function, <:AbstractVector{<:Symbol}}, 
#										 digits::ODict, 
#										 allparams::UODict)
#
#		@assert NrParamSets==1 
#
#		return new(NrParamSets,
#							 typical_allparams(uk, allparams),
#							 construct_get_fname(path, uk, digits)
#							 )
#
#	end 
#
#	function ParamFlow(path, 
#										 uk::Union{<:Function, <:AbstractVector{<:Symbol}}, 
#										 digits::ODict, 
#										 allparams::UODict)
#
#		ParamFlow(1, path, uk, digits, allparams)
#
#	end 

#	function ParamFlow(path, tuples::Vararg{<:Tuple})
#
#
#	end  
#
#
#
#
#	function ParamFlow(path, tuples::Vararg{<:Tuple})
#
#
#	end
#
#	function ParamFlow(N::Union{<:Nothing,<:Int}, allparams::Function, 
#										 FNG::FilenameGenerator)::ParamFlow
#
#		new(N, allparams, (getproperty(FNG, p) for p in propertynames(FNG))...)
#	
#	end 
#
#
#
#	function ParamFlow(allparams::Function, FNG::FilenameGenerator)::ParamFlow
#
#		ParamFlow(nothing, allparams, FNG)
#
#	end 
#
#
#
#	function ParamFlow(allparams::Function, arg1, arg2, args...)::ParamFlow
#	
#		ParamFlow(allparams, FilenameGenerator(arg1, arg2, args...))
#	 
#	end
#
#
#	function ParamFlow(N::Union{<:Nothing, <:Int}, allparams::Function, 
#										 arg1, arg2, args...)::ParamFlow
#
#		ParamFlow(N, allparams, FilenameGenerator(arg1, arg2, args...))
#
#	end 
#
#
#
#
#	function ParamFlow(NrParamSets::Int, input_dict::UODict, args...)::ParamFlow
#	
#		FNG = FilenameGenerator(args...) 
#	
#		function allparams(args::Vararg{<:UODict,N}) where N 
#	
#			N < NrParamSets || error()  
#
#			return typical_allparams(FNG.usedkeys, input_dict)
#	
#		end 
#
#		return ParamFlow(NrParamSets, allparams, FNG)
#	 
#	end
#


#	function ParamFlow(M::Module, input_dict::UODict, args...)::ParamFlow
#
#		ParamFlow(M.NrParamSets, input_dict, args...)
#
#	end 
#
#
#	function ParamFlow(input_dict::UODict, args...)::ParamFlow 
#
#		ParamFlow(1, input_dict, args...)
#
#	end 
#
#




#	function ParamFlow(N::Union{<:Nothing,<:Int},
#										 allparams::Function, arg1, arg2, args...)::ParamFlow
#
#		ParamFlow(N, allparams, FilenameGenerator(arg1, arg2, args...))
#
#	end 

end 


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#




#
#function get_usedkeys(S::Union{<:FilenameGenerator, <:ParamFlow}, args...)
#
#	get_usedkeys(S.usedkeys, args...)
#
#end 
#
#
#function get_params_digits(S::Union{<:FilenameGenerator, <:ParamFlow}, 
#													 args...)
#
#	get_params_digits(S.params_digits, args...)
#
#end 
#
#


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function good_methods(NrParamSets::Int, F::Function, args...)::Bool

	possib_types = [UODict, 
									AbstractVector{<:UODict}, 
									NTuple{N,<:UODict} where N]

	any(Base.product((possib_types for _=1:NrParamSets)...)) do tup
	
		hasmethod(F, tup, args...)

	end || return false 


	return true 

end 

good_methods(::Nothing, args...)::Bool = good_methods(1, args...)

good_methods(PF::ParamFlow, a...)::Bool = good_methods(PF.NrParamSets, a...) 




#===========================================================================#
#
# adds kwargs, if F accepts them: those in 'new_kwargs' and 'kwargs'
#
#---------------------------------------------------------------------------#


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

	return Utils.add_args_kwargs(F; kwa1..., kwa2...)

end 




#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#



struct Calculation

	name::String

	PF::ParamFlow
	
	Compute::Function 

	Read::Function 

	FoundFiles::Function 


	function Calculation(name::Union{<:AbstractString,<:Symbol}, PF::ParamFlow,
											 CRF::Vararg{<:Function,3}; kwargs...)

		new(string(name), 
				PF, (add_kwargs(PF, f; kwargs...) for f in CRF)...)

	end 


	function Calculation(name::Union{<:AbstractString,<:Symbol}, PF::ParamFlow,
											 Compute::Function; kwargs...)

		C = add_kwargs(PF, Compute; kwargs...) 

		false_(args...;kwargs...)::Bool=false 

		return new(string(name), PF, C, C, false_)

	end 



	function Calculation(name::Union{<:AbstractString,<:Symbol}, PF::ParamFlow, 
											 Compute::Function, ::Nothing, ::Any; kwargs...)

		Calculation(name, PF, Compute; kwargs...) 

	end 


	function Calculation(name::Union{<:AbstractString,<:Symbol}, PF::ParamFlow, 
											 Compute::Function, ::Nothing, ::Nothing; kwargs...)

		Calculation(name, PF, Compute; kwargs...)

	end 


	function Calculation(PF::ParamFlow, Compute::Function, args...; kwargs...)

		Calculation("", PF, Compute, args...; kwargs...)

	end 




	function Calculation(name::Union{<:AbstractString,<:Symbol}, 
											 PF::ParamFlow, M::Module; kwargs...)

		Calculation(name, PF, Utils.getprop(M, [:Compute,:Read,:FoundFiles])...;
								kwargs...)

	end 


	function Calculation(PF::ParamFlow, M::Module; kwargs...)

		Calculation(tostr(M), PF, M; kwargs...)

	end 



end 



tostr(C::Calculation)::String = C.name  

get_NrPSets(m::Module)::Int = m.NrParamSets
get_NrPSets(pf::ParamFlow)::Int = pf.NrParamSets 
get_NrPSets(c::Calculation)::Int = get_NrPSets(c.PF)


get_allP(m::Module, args...) = m.allparams(args...) 
get_allP(pf::ParamFlow, args...) = pf.allparams(args...) 
get_allP(m::Calculation, args...) = get_allP(m.PF, args...)



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


#get_allP, get_NrPSets -- applied to M

				# FIRST CALL 

function convert_params(M::Union{<:Module, <:ParamFlow, <:Calculation};
												repl=nothing,
												kwargs...)
	convert_params(M, [];
								 repl=combine_functions_addrem(repl, get_PF_kwargs(M, :adjust, nothing)),
								 kwargs...)

end  



					# FURTHER CALLS 
function convert_params(M::Union{<:Module, <:ParamFlow, <:Calculation},
												args::AbstractVector; 
												repl::Function, 
												get_new=(l,a)->get_allP(M, a...),
												convert_one=identity,
												convert_many=map,
												final_f=identity)

	for level in 1:get_NrPSets(M)

		push!(args, repl(level, args..., get_new(level, args)))

	end


	return map(args) do P

		if applicable(convert_one, P)
			
			return convert_one(P)

		elseif Utils.isList(P) && all(p->applicable(convert_one, p), P)

			return convert_many(convert_one, P)

		else 

			error("\n $P \nParameter structure not understood.")

		end 

	end |> final_f 


end

#===========================================================================#
#
# Particularizes convert_params in order to obtain a single dictionary
# 		with unique, string keys, from a structured set P
#
# 		(the inverse operation of convertParams_fromPlot)
#
#---------------------------------------------------------------------------#

function convertParams_toPlot(M::Union{<:Module, <:ParamFlow, <:Calculation},
															P::Utils.List; kwargs...)

	isempty(P) && return convertParams_toPlot(M; kwargs...)

	return convertParams_toPlot(M; get_new = (level,args)->P[level], kwargs...)

end 


function convertParams_toPlot(M::Union{<:Module, <:ParamFlow, <:Calculation},
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

function convertParams_fromPlot(M::Union{<:Module,<:ParamFlow,<:Calculation},
																P::UODict; kwargs...)

	function repl(params::UODict)::OrderedDict

		K = Utils.mapif(k1->(k1,convertKey_toPlot(params,k1)),
		
										k -> haskey(P,k[2]),
		
										collect(keys(params)))
		
		return OrderedDict(k1 => P[k2] for (k1,k2) in K)

	end 

	repl(P::Utils.List)::Vector{OrderedDict} = repl.(P)

	repl(l::Int, P...) = repl(P[l])


	return convert_params(M; repl=repl, kwargs...)

end


##################### REFORMULTE WITH Utils.Backtracking !! 

#===========================================================================#
#
#	Fill internal parameters
#
#---------------------------------------------------------------------------#


function fillParams_internalP(M::Union{<:Module, <:ParamFlow, <:Calculation},
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



function fillParams_internalPs(M::Union{<:Module,<:ParamFlow,<:Calculation}, 
															 P::Utils.List, ::Nothing;
															kwargs...)
	P

end 


function fillParams_internalPs(M::Union{<:Module,<:ParamFlow,<:Calculation}, 
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

get_PF_kwargs(M::Module, args...) = Utils.getprop(M, args...)
get_PF_kwargs(pf::ParamFlow, args...) = Utils.getprop(pf, args...)
get_PF_kwargs(c::Calculation, args...) = get_PF_kwargs(c.PF, args...)


function get_paramcombs(M::Union{<:Module, <:ParamFlow, <:Calculation};
												cond=nothing, repl=nothing) # First call

	get_paramcombs(M, 1, [];
								 cond = combine_functions_cond(cond, 
																get_PF_kwargs(M, :constraint, nothing)),
								 repl = combine_functions_addrem(repl, 
																get_PF_kwargs(M, :adjust, nothing)),
								)
end 



function get_paramcombs(M::Union{<:Module, <:ParamFlow, <:Calculation},
												level::Int, old::AbstractVector; 
												cond::Function, repl::Function) # Further calls 

	#M must have M.allparams::Function, M.NrParamSets::Int,

	news = Utils.mapif(pcomb -> vcat(old, [repl(level, old...,  pcomb)]),

										 new_ -> cond(level, new_...),

										 Utils.AllValCombs(repl(level, old..., get_allP(M, old...))))




	level == get_NrPSets(M) && return news

	return Utils.flatmap(news) do new 
				
						get_paramcombs(M, level+1, new; cond=cond, repl=repl) 
	end

end




function f_get_plotparams(M::Union{<:Module, <:Calculation}, 
													rmv_internal_key::Nothing=nothing)

#	getpp() = convertParams_toPlot(M) 

	getpp(P...) = convertParams_toPlot(M, P)
												
#	return getpp 

end 


function f_get_plotparams(M::Union{<:Module, <:Calculation}, 
													rmv_internal_key::Function)

	getpp(P...) = convertParams_toPlot(M, P; repl=rmv_internal_key) 

end 



function f_fill_internal(M::Union{<:Module, <:Calculation}, args...)
	
	fill_int(P::Utils.List) = Parameters.fillParams_internalP(M, P, args...)
	
end 


#############################################################################
end
