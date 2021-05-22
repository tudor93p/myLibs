module Parameters
#############################################################################

import ...OrderedDict, ...Utils 


const ODict = Union{<:OrderedDict, <:NamedTuple}
const UDict = AbstractDict
const UODict = Union{<:AbstractDict, <:OrderedDict, <:NamedTuple}










#############################################################################
module Operations
#############################################################################

import ..UODict, ...Utils

import Random 



#===========================================================================#
#
# applies several add or remove functions recursively
#
#---------------------------------------------------------------------------#


function combine_functions_addrem(args...)::Function

	good_fs = Utils.flat(args...; keep=x->isa(x,Function))


	function combined_addrem(l::Int, p::T)::T where T<:UODict

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

	return function combined_cond(P::T, l::Int)::Bool where T<:UODict

		for cond in good_fs
			
			!cond(P,l) && return false 

		end 

		return true 

	end


end 





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

function convert_params(M::Module, args=[];
												get_new=(l,a)->M.allparams(a...),
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


	for level in 1:M.NrParamSets

		new = get_new(level,args)

		if isnothing(repl) push!(args,new) else 

				push!(args, act_on_plev(a->repl(level,a))(new))

		end

	end


	return final_f(map(act_on_plev(convert_one),args))

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



#===========================================================================#
#
# Particularizes convert_params in order to obtain a single dictionary
# 		with unique, string keys, from a structured set P
#
# 		(the inverse operation of convertParams_fromPlot)
#
#---------------------------------------------------------------------------#

function convertParams_toPlot(M::Module, P=nothing; kwargs...)

	!isnothing(P) && return convertParams_toPlot(M;
																	 get_new = (level,args)->P[level],
																	 kwargs...)

	convert_params(M;

		convert_one = function (params)

										#map(collect(params)) do (k,v)

										#	convertKey_toPlot(params,k) => v

										#end

										[convertKey_toPlot(params,k)=>v for (k,v) in params]


									end,

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

function convertParams_fromPlot(M::Module,P;kwargs...)

	convert_params(M;

		repl = function (level,params)

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


function fillParams_internalP(M::Module, P, add_internal_param;
															kwargs...)

	isnothing(add_internal_param) && return P

	return convert_params(M;
	
						get_new = (level,args) -> P[level],
	
						repl = add_internal_param,

						kwargs...
	
					)
end 



function fillParams_internalPs(M::Module, P, add_internal_params;
															kwargs...)

	isnothing(add_internal_params) && return P

	return map(add_internal_params) do add

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

function get_paramcombs(M::Module, level=1, old=[];
												cond=nothing, repl=nothing)

	repl_ = combine_functions_addrem(repl)
	cond_ = combine_functions_cond(cond)

	raw = repl_(level,M.allparams(old...))
	
	news = Utils.mapif(pcomb -> vcat(old...,pcomb),

										 isnothing(cond) || new->cond_(new, level),

										 Utils.AllValCombs(raw))

	level == M.NrParamSets && return news

	return Utils.flatmap(news) do new 
				
						get_paramcombs(M, level+1, new; cond=cond, repl=repl_) 
	end

end




function f_get_plotparams(M, rmv_internal_key=nothing)

	getpp(P=nothing) = convertParams_toPlot(M, P; repl=rmv_internal_key)

end 


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

function replace_parameter_fs(getP, Keys, level=1)::Tuple{Function,Function}

	isa(level,Int) || Utils.isList(level, Int) || error()

	islev = in(Utils.flat(level))

	K = Utils.flat(Keys)

	rmv(l,p) = islev(l) ? filter(kv->!in(kv.first, K), p) : p

	if getP isa Function 
		
		add(l,p) = islev(l) ? merge(p, Dict(zip(K, getP(l,p)))) : p

		return rmv,add

	end

	return rmv,(l,p)->p
	#return Dict(:rmv_internal_key=>rmv, :add_internal_param=>add)

end 


#############################################################################
end
#############################################################################






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

	usedkeys::Union{<:Nothing, <:Vector{<:Symbol}, <:Function}
	
	params_digits::Function 

	path::String  

	get_fname::Function

end 


function get1(M::T, ks::Vararg{<:Symbol}) where T 

	T<:Module || return M

	i = findfirst(in(names(M, all=true)), ks)
							
	!isnothing(i) && return getproperty(M, ks[i]) 

	error("None of $ks found in module $M")

end 


# ----- # 

function FilenameGenerator(usedkeys::Union{<:Nothing, 
																	<:Vector{<:Symbol}, 
																	<:Function}, 
									params_digits::Function, 
									path::String)::FilenameGenerator

	function get_fname(args::Vararg{<:UODict})::Function
		
		fname(tostr([path, tostr(params_digits, args...)]))

	end  

	return FilenameGenerator(usedkeys, params_digits, path, get_fname)

end 





function FilenameGenerator(usedkeys, params_digits::Function, path)::FilenameGenerator

	FilenameGenerator(get1(usedkeys, :usedkeys), params_digits, 
					 prefix(get1(path,:path)))

end 

function FilenameGenerator(usedkeys, args_params_digits::Tuple, path)::FilenameGenerator

	FilenameGenerator(get1(usedkeys, :usedkeys), 
					 typical_params_digits(args_params_digits...), path)

end 

function FilenameGenerator(usedkeys::Union{<:AbstractVector{<:Symbol}, <:Function},
									digits::ODict, path)::FilenameGenerator

	FilenameGenerator(usedkeys, (usedkeys, digits), path)

end  

function FilenameGenerator(M_usedkeys::Module, digits::ODict, path)::FilenameGenerator

	FilenameGenerator(get1(M_usedkeys, :usedkeys), digits, path)
	
end  


function FilenameGenerator(usedkeys, M_digits::Module, args...)::FilenameGenerator

	FilenameGenerator(usedkeys, get1(M_digits, :digits, :params_digits), args...)

end 
#
#	function FilenameGenerator(usedkeys, params_digits, M_path::Module)::FilenameGenerator
#
#		FilenameGenerator(usedkeys, params_digits, get1(M_path, :path))
#
#	end 








function FilenameGenerator(params_digits::Union{<:Tuple, <:Function}, path)::FilenameGenerator

	FilenameGenerator(nothing, params_digits, path)

end 

function FilenameGenerator(M_params_digits::Module, path)::FilenameGenerator

	FilenameGenerator(get1(M_params_digits, :params_digits), path) 

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

tostr(M::Module)::String = tostr(Base.fullname(M)[2:end]) # exclude :Main




#function prefix(path::Union{<:AbstractString,<:Utils.List}, func::Function, 
#							 kwargs...)::String where N 
#
#	prefix(path; kwargs...)
#
#end 


function tostr(func::Function, arg::UODict)::String

	params_to_string(func(arg))

end 


function tostr(func::Function, args::Vararg{<:UODict}; kwargs...)::String

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

#function typical_params_digits(M::Module, args...)::Function
#
#	typical_params_digits(M.usedkeys, args...)
#
#end



function params_digits_(P::Tp, usedkeys::Function, digits::Td
											 )::Tuple{Tp, AbstractVector{Symbol}, Td} where {
																								Tp<:UODict, Td<:ODict}

	params_digits_(P, usedkeys(P), digits)

end 


function params_digits_(P::Tp, usedkeys::Tu, digits::Td
											)::Tuple{Tp, Tu, Td} where {Tp<:UODict,
																									Tu<:AbstractVector{Symbol},
																									Td<:ODict}

	(P, usedkeys, digits)

end 




function typical_params_digits(usedkeys::Union{<:AbstractVector{Symbol},
																							 <:Function}=Symbol[],
															 digits::ODict=OrderedDict())::Function

	(P::UODict) -> params_digits_(P, usedkeys, digits)

end


#function merge_params_digits(args::Vararg{Any,N}) where N
#	
#	N==1 && return merge_params_digits(args...)
#				# apparently one argument is given, must be in fact a group
#				
#	function result(P...) 
#
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
#		out(tuples) = Tuple(map(enumerate([merge, union, merge])) do (i,f)
#	
#													f([t[i] for t in tuples]...)
#
#												end)
#
#
##		tuples_ = ((p1,u1,d1), (p2,u2,d2) etc)
#
#		tuples_[1][1] isa Union{AbstractDict,NamedTuple} && return out(tuples_)
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
#
#
#
#





#############################################################################
#end
#############################################################################

struct ParamFlow

	allparams::Function

	usedkeys::Union{<:Nothing, <:Vector{<:Symbol}, <:Function}
	
	params_digits::Function 

	path::String  

	get_fname::Function

	function ParamFlow(allparams::Function, FNG::FilenameGenerator)::ParamFlow
	
		new(allparams, (getproperty(FNG, p) for p in propertynames(FNG))...)
	
	end 
	
	function ParamFlow(allparams::Function, arg1, arg2, args...)::ParamFlow
	
		ParamFlow(allparams, FilenameGenerator(arg1, arg2, args...))
	 
	end
	
	function ParamFlow(input_dict::UODict, args...)::ParamFlow 

		ParamFlow((1,input_dict), args...)

	end 
	
	function ParamFlow((NrParamSets,input_dict)::Tuple{<:Int,<:UODict}, 
										 args...)::ParamFlow
	
		FNG = FilenameGenerator(args...) 
	
		function allparams(args::Vararg{<:UODict,N}) where N 
	
			N < NrParamSets || error() 
			return Operations.typical_allparams(FNG.usedkeys, 
																					input_dict[:allparams])
	
		end 
	
		return ParamFlow(allparams, FNG)
	 
	end
	

end 

#usedkeys = ...
#NrParamSets = ...
#
#DATAROOT defined in Constants 
#
#function get_PF(input_dict::AbstractDict)::FilenameGenerator

#ParamFlow((NrParamSets, input_dict), usedkeys, input_dict, DATAROOT)
#end 




#############################################################################
end
