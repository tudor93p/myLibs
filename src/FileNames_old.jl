module FileNames
#############################################################################

import ..Utils

const ROOT = "Data"


using OrderedCollections: OrderedDict
#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

# fname("sometask",(libA,pA1,pA2),(libB,pB),(libC,pC1,pC2,pC3,pC4))
# path = "root/sometask/signature(pA1)/signature(pA2)/sinature(pB)/..."






function fname(name::String, func::Function, arg; ignore=[])::Function

	strings = Utils.flatmapif(params_to_string, !isempty, [func(args)])
	
	prefix = mkpath(join([ROOT;name; strings],"/"))

	return (x="") -> prefix*"/"*x

end



function fname(name::String, func::Function, args...; ignore=[])::Function

	strings = Utils.flatmapif(params_to_string, !isempty, func(args...))


	prefix = mkpath(join([ROOT; name; strings],"/"))

	return (x="") -> prefix*"/"*x

end



function fname(name::Union{AbstractVector,Tuple}, args...; 
							 																		kwargs...)::Function
	fname(join(name,"/"), args...; kwargs...)

end

function fname(name::Union{AbstractVector,Tuple}; kwargs...)::Function

	(args...; kw1...)-> fname(join(name,"/"), args...; kwargs..., kw1...)

end


#function fname(name::AbstractVector, args...; kwargs...)

#	fname(join(name,"/"), args...; kwargs...)

#end


function fname(M::Module, n::Int64=1,args...; kwargs...)::Function

	fname(M, n, M.params_digits, args...; kwargs...)

end


function fname(M::Module, n::Int64=1, params_digits::Function, args...; kwargs...)::Function

	fname(Base.fullname(M)[end-n+1:end], params_digits, args...; kwargs...)

end






#===========================================================================#
#
# Turn one set of parameters into a string
#
#---------------------------------------------------------------------------#


function params_to_string(params,usedkeys,Digits;separator='-')

	K = intersect(keys(Digits),intersect(keys(params),usedkeys))

	isempty(K) && return ""

  return join([Utils.nr2string(params[k],Digits[k]) for k in K],separator)

end

params_to_string((p,u,d);kwargs...) = params_to_string(p,u,d;kwargs...)







#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function typical_params_digits(M::Module, args...)::Function

	typical_params_digits(M.usedkeys, args...)

end

function typical_params_digits(usedkeys::AbstractVector{Symbol}=[], 
															digits::T=Dict())::Function where T<:Union{NamedTuple,OrderedDict}

	(P::AbstractDict) -> (P, usedkeys, digits)

end


function merge_params_digits(args::Vararg{Any,N}) where N
	
	N==1 && return merge_params_digits(args...)
				# apparently one argument is given, must be in fact a group
				
	function result(P...) 


		tuples_ = map(args) do arg 
			
			typeof(arg)<:Module && return arg.params_digits(P...)
			
			typeof(arg)<:Function && return arg(P...)
			
			return arg
			
		end


		out(tuples) = Tuple(map(enumerate([merge, union, merge])) do (i,f)
	
													f([t[i] for t in tuples]...)

												end)


#		tuples_ = ((p1,u1,d1), (p2,u2,d2) etc)

		tuples_[1][1] isa Union{AbstractDict,NamedTuple} && return out(tuples_)

	
#		tuples_ = ( [ (pi1,ui1,di1) etc ], [ (pj1,uj2,dj2) etc ], etc  )

		return [out(tuples) for tuples in zip(tuples_...)]

#		out( (pi1,ui1,di1), (pj1,uj1,dj1) ), out( (pi2,ui2,...), (j2) )

	end



	for arg in args

		typeof(arg)<:Module && continue

		typeof(arg)<:Function && continue

		(p, u, d) = arg
		
		return result(p)

#		if at least one tuple (p, u ,d) is given, will return the total tuple 
#								(total_p, total_u, total_d)
			
	end

	return (P...) -> result(P...)
	

	
# if only modules are given, 
#						 return the function which computes the total tuple

end










#############################################################################
end

