module GFD
#############################################################################

import ..MetaDiGraph




import ..Graph
import ..LayeredSystem: islayer, islead,
												node_right, node_left, 
												get_lead_GFf, get_graphH,
												get_node


import ..GF,..SelfEn

#===========================================================================#
#
#	Helper functions for dealing with the physical system stored in the graph
#
#---------------------------------------------------------------------------#


function next_node(g::MetaDiGraph, 
									 n::Int64, 
									 d::AbstractString)::Union{Nothing,Int}

	d=="right" && return node_right(g, n)

	d=="left" && return node_left(g, n)

	error("d should be left or right")

end 

meets_layer(::MetaDiGraph, ::Nothing, ::AbstractString="")::Bool = false  

meets_layer(::MetaDiGraph, ::Int)::Bool = true  

function meets_layer(g::MetaDiGraph, n::Int, d::AbstractString)::Bool

	(d=="both" || islayer(g,n)) && return true

	return meets_layer(g, next_node(g,n,d), d)

end

function dir_layer(g::MetaDiGraph, n::Int)::String 

	for d in ("right","left")

		meets_layer(g,n,d) && return d 
		
	end 
	
	error()

end 

function dir_out(g::MetaDiGraph, n::Int)::String 

	for d in ("right","left")

		meets_layer(g,n,d) || return d 
		
	end 
	
	error()

end 

dir_layer((g,)::Tuple, n::Int)::String = dir_layer(g, n)


# ---- coupling node(s) left or right ----- #
	


function bring_closer(g::MetaDiGraph, n::Int, m::Int
											)::Tuple{Bool,Int,Int}

	if isempty(Graph.FindPath(g,n,m))  #means m<n
		
		return (false, node_left(g,n), node_right(g,m))

	else 

		return (true, node_right(g,n), node_left(g,m))

	end

end  



Gnm_inds((g,)::Tuple, args...) = Gnm_inds(g, args...)

function Gnm_inds(g::MetaDiGraph, n::Int, m::Int, d::AbstractString
						 )::NTuple{2,Tuple{Int,Int,String}}

	nisleft, n1, m1  = bring_closer(g, n, m)
	
	if d in ("both","left")

		return nisleft ? ((n,m1,"left"), (m,m,d)) : ((n,n,d), (n1,m,"left"))  
	elseif d=="right"

		return nisleft ? ((n,n,d), (n1,m,d)) : ((n,m1,d), (m,m,d)) 

	end 

end  

function pure_lead((g,)::Tuple, n::Int, m::Int, d::AbstractString)::Bool 
	
	n==m && islead(g,n) && !meets_layer(g,n,d)

end 

function cpl_lead((g,)::Tuple, n::Int, dir::AbstractString)::Bool
	
	islead(g,n) && dir=="both"

end 



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function H((g,)::Tuple{MetaDiGraph, Number, Vararg}, args...
					)::AbstractMatrix{ComplexF64}

	get_graphH(g, args...)

end 

function H((g,data_H,)::Tuple{MetaDiGraph, Dict, Vararg}, args...
					)::AbstractMatrix{ComplexF64}

	get_graphH(g, data_H, args...)

end  

function Energy((_, _, E,)::Tuple{MetaDiGraph, Dict, T, Vararg}
							 )::T where T<:Number  
	E 
end  

function Energy((_, E,)::Tuple{MetaDiGraph, T, Vararg})::T where T<:Number 
	E
end 

function reEnergy((_, _, E, rE)::Tuple{MetaDiGraph, Dict, Number, T, Vararg}
								 )::T where T<:Number  
	rE 

end 


function reEnergy((_, E, rE)::Tuple{MetaDiGraph, T, Vararg}
								 )::T where T<:Number 

	rE 

end 

function coupling_inds(g::MetaDiGraph, 
											 src::Int, 
											 dir::AbstractString
														 )::Vector{Tuple{Int,Int,String}}

	if dir in ("right","left")

		 dst = next_node(g, src, dir) 

		 return isnothing(dst) ? Tuple{Int,Int,String}[] : [(src,dst,dir)]

	end 

	return vcat((coupling_inds(g, src, d) for d=("right","left"))...)

end  



function coupling_inds((g,)::Tuple, args...)::Vector{Tuple{Int,Int,String}}

	coupling_inds(g, args...)
 
end  





#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function lead_GF!((_,storage,)::Tuple, data::Tuple, args... 
									)::AbstractMatrix{ComplexF64}

	lead_GF!(storage, data, Energy(data), args...)

end  

function lead_GF!(storage::AbstractDict,
											(g,)::Tuple{MetaDiGraph,Vararg},
											E::Number,

								lead_label::AbstractString,
								lead_uc::Int,
								args...
							 )::Matrix{ComplexF64} 

	out = get!(storage, lead_label) do 

		get_lead_GFf(g, lead_label)(E)

	end  

	return out[min(lead_uc,end)]

end  

function lead_GF!(storage::AbstractDict,
											data::Tuple,
											E::Number,
								n::Int,
								args...
							 )::Matrix{ComplexF64} 

	lead_GF!(storage, data, E, Graph.get_prop(data[1], n, :name)...)

end  



function sys_GF!(storage::Tuple,
						data::Tuple, 
						n::Int, m::Int, d::AbstractString 
					)::AbstractMatrix{ComplexF64}

	get!(storage[1], (n,m,d)) do 
											Gnm!(storage, data, n,m,d) end 

end  



#===========================================================================#
#
# 	Green's function using decimation method for at most two leads.
#
#		Implementation as explained by Lewenkopf+Mucciolo in
#
#		"The recursive Green’s function method for graphene"
#
#				J Comput Electron (2013) 12: 203–231
#
#---------------------------------------------------------------------------#

function coupling!(storage::Tuple, data::Tuple, args...)::Base.Generator

	((H(data,m,n),
		G!(storage,data,m,m,d)) for (n,m,d)=coupling_inds(data,args...))

end 
		
function GF_EmH!(storage::Tuple, data::Tuple, n::Int, d::AbstractString
								)::AbstractMatrix{ComplexF64} 

	GF(reEnergy(data), H(data, n), coupling!(storage, data, n, d)...) 

end   

function GF_gS!(storage::Tuple,
								data::Tuple,
								n::Int, d::AbstractString
								)::AbstractMatrix{ComplexF64} 

	GF(lead_GF!(storage, data, n),
		 coupling!(storage, data, n, dir_layer(data, n))...)

end 




function G!(storage::Tuple, args...)::AbstractMatrix{ComplexF64}

	if pure_lead(args...)
		lead_GF!(storage, args...)
	else 
		sys_GF!(storage, args...)
	end 

end

function Gnm!(storage::Tuple, data::Tuple,  
						 n::Int, m::Int, d::AbstractString,
						 )::AbstractMatrix{ComplexF64}

	n==m && return (cpl_lead(data,n,d) ? GF_gS! : GF_EmH!)(storage, data, n, d)

	(a1,b1,d1),(a2,b2,d2) = Gnm_inds(data, n, m, d) 

	return *(G!(storage, data, a1,b1,d1), 
					 H(data, b1, a2), 
					 G!(storage, data, a2,b2,d2))

end 


#===========================================================================#
#
# User methods 
#
#---------------------------------------------------------------------------#


function G!(storage::Tuple, data::Tuple, 
						(n1,i1,n2,i2)::Tuple{String,Int,String,Int},
						dir::AbstractString,
						::Nothing=nothing,
						 )::AbstractMatrix{ComplexF64}

	G!(storage, data, get_node(data[1], n1, i1), get_node(data[1], n2, i2), dir)

end 

function G!(storage::Tuple, data::Tuple, 
						ns::Tuple, dir::AbstractString, slicer::Function,
						 )::AbstractMatrix{ComplexF64}

	Ns, inds = slicer(ns...)

	return view(G!(storage, data, Ns, dir), inds...)

end 



function unpack_argsG(ni1ni2::Tuple{Tuple,Tuple}; kwargs...
							)::Tuple{Tuple{String,Int,String,Int},String}

	unpack_argsG(ni1ni2...; kwargs...)

end  

function unpack_argsG(n1i1::Tuple{<:AbstractString,Int},
											n2i2::Tuple{<:AbstractString,Int}=n1i1;
											kwargs...
							)::Tuple{Tuple{String,Int,String,Int},String}

	unpack_argsG(n1i1..., n2i2...; kwargs...)

end 


function unpack_argsG(name1::AbstractString, index1::Int64,
						 	name2::AbstractString=name1, index2::Int64=index1;
							dir::AbstractString="both"
							)::Tuple{Tuple{String,Int,String,Int},String}

	(name1,index1,name2,index2), dir 

end  





function get_G!(storage::NTuple{2,Dict}, 
								data::Tuple{MetaDiGraph, Vararg},
								slicer::Union{Function,Nothing}=nothing,
								)::Function

	function _G!(a...; k...)::AbstractMatrix{ComplexF64}

		G!(storage, data, unpack_argsG(a...; k...)..., slicer)

	end 

end





#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#




function get_SE!(storage::NTuple{2,Dict},
								 data::Tuple{MetaDiGraph, Vararg},
								 slicer::Function,
								)::Function


	function self_en(k::AbstractString, uc::Int=2)::Matrix
	
		(K,),(s,) = slicer(k,uc) 
		
		i = max(uc,2) 
		
		return SelfEn(view(H(data, K, i, K, i-1), s, s),
									view(lead_GF!(storage, data, K, i), s, s))

	end 

end 	









































































#############################################################################
end  # module GFD 
#############################################################################

