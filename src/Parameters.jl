module Parameters
#############################################################################

import ..OrderedDict, ..Random 

import ..Utils




#===========================================================================#
#
# applies several add or remove functions recursively
#
#---------------------------------------------------------------------------#


function combine_functions_addrem(args...)

	good_fs = Utils.flat(args...; keep=x->isa(x,Function))

	function combined_addrem(l::Int, p::T) where T
		
		T<:Union{AbstractDict,NamedTuple} || return combined_addrem.(l,p)

		p1 = copy(p)

		for fi in good_fs

			p1 = fi(l, p1)

		end 

		return p1

	end

	return combined_addrem

end 


function combine_functions_cond(args...)

	good_fs = Utils.flat(args...; keep=x->isa(x,Function))

	return function combined_cond(P, l::Int)

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

function convertKey_toPlot(p,k;separator="__")

	:Label in keys(p) && return join([p[:Label], k], separator) 

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


function typical_allparams(usedkeys, allp)

	filter(p->p.first in usedkeys, OrderedDict{Symbol,Any}(pairs(allp)))

end


#function merge_allparams(Ms::Vararg{Module,N}; args=[], kwargs...) where N
function merge_allparams(inputs...; args=[], kwargs...)


	out = map(inputs) do i
	
						typeof(i)<:Module && return i.allparams(args...; kwargs...)
				
						typeof(i)<:AbstractDict && return i
																# already the result of typical_allparams
			
						return typical_allparams(i...)
				
					end
					
#	Utils.isList(out,AbstractDict)

	eltype(out) <: Union{NamedTuple,AbstractDict} && return merge(out...)

	if Utils.isList(eltype(out),Union{NamedTuple,AbstractDict})

		return [merge(o...) for o in zip(out...)]

	end

	println()

	println(out)


	error("Element type not undestood")

end





#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function replace_parameter_fs(getP, Keys, level=1)

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
