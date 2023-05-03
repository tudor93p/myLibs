module MeshInPlace
############################################################################

import Distributed: @spawnat, myid, workers, nworkers, procs
#import DistributedArrays: DArray, SubOrDArray, localindices, localpart, dzeros 

using SharedArrays:SharedArray


import ..Utils 


const SubOrArray{T,N} = Union{Array{T,N}, SubArray{T,N,<:Array}}
const SubOrSArray{T,N} = Union{SharedArray{T,N}, SubArray{T,N,<:SharedArray}}

const ENDPOINT::Bool = true 

#===========================================================================#
#
# Small tools 
#
#---------------------------------------------------------------------------#


function check_parallel_possible(;parallel::Bool=false, kwargs...)::Bool

	parallel && nworkers()>1 

end  


"""
Obtain args and kwargs from `akfun` and call `fun`
"""
function wrapper_ak(akfun::Function, fun::Function, args...; kwargs...) 

	a, k = akfun(args...; kwargs...)

	return fun(a...; k...) 

end  

function wrapper_ak(akfun::Function, 
										funs::Tuple{Vararg{<:Function}},
										args...; kwargs...) 

	a, k = akfun(args...; kwargs...)

	return [f(a...; k...) for f in funs]

end 


#function wrapper_scalar_inplace(f::Function)::Function 
#
#	f!(V::AbstractVector{<:Number}, args...; kwargs...) 
#
#
#f(args...; kwargs...)
#









#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#



""" 
Get mesh dimension `dim` along direction `dir` 

To be used in `selectdim(storage, dim, i)`, `axes(storage, dim)` etc.,
compatible with [`_storage_size`](@ref).

# Arguments 
- `Nm::Int`: dimensionality of the mesh 
- `Ns::Int`: `ndims(storage)>=Nm` 
"""
function get_mesh_dim(Ns::Int, Nm::Int, dir::Int)::Int 

	@assert Ns>=Nm && 1<=dir<=Nm 

	Ns-Nm+dir 

end 

function get_mesh_dim(::Union{NTuple{Ns,<:Union{Int,UnitRange{Int}}},
															 AbstractArray{<:Number,Ns}
															 },
											args::Int...
											)::Int where Ns 

	get_mesh_dim(Ns, args...)

end  




function mesh_axes(ns::Vararg{Int,Nm}; kwargs...)::NTuple{Nm,Int} where Nm

	mesh_axes(ns; kwargs...) 

end 

function mesh_size(ns::Vararg{Int,Nm}; kwargs...)::NTuple{Nm,Int} where Nm

	mesh_size(ns; kwargs...) 

end 

function mesh_size(ns::NTuple{Nm,Int}, args...;
									 endpoint::Bool=ENDPOINT,
									 kwargs...)::NTuple{Nm,Int} where Nm

	endpoint ? ns : ns .- 1 

end  

function mesh_axes(ns::NTuple{Nm,Int}, args...; 
									 kwargs...)::NTuple{Nm,Base.OneTo} where Nm

	map(Base.OneTo, mesh_size(ns; kwargs...))

end 

"""
Size of a storage array. 

# Arguments
- `s::NTuple{Na,Int}`: size of the array `A` at each mesh point, `ndims(A)==Na` 
- `ns::NTuple{Nm,Int}`: mesh size, dimensionality `Nm`

# Keyword arguments 
- `endpoint::Bool=true`: whether the last and first mesh points differ
"""
function _storage_size(
											 s::Tuple{Vararg{Int}},
											 ns::Tuple{Vararg{Int}},
											args...;
											kwargs...
											)::Tuple{Vararg{Int}}

	Utils.tuplejoin(isempty(s) ? (1,) : s, mesh_size(ns; kwargs...))

end 





#===========================================================================#
#
# _prep_init_array and _init_storage take identical same args and kwargs
# internal methods 
#
#---------------------------------------------------------------------------#


function _prep_init_array(
											s::Tuple{Vararg{Int}},
											ns::NTuple{Nn,Int},
											args...;
											distr_axis::Int=Nn,
											kwargs...
											)::Tuple{Tuple{Vararg{Int}}, Vector{Int}, Vector{Int}
															 } where Nn

	array_size = _storage_size(s, ns; kwargs...)

	Na = length(array_size) - Nn

	d = ones(Int, length(array_size))


	parallel = check_parallel_possible(;kwargs...)


	parallel || return (array_size, [myid()], d) # just the master proc

	@assert 1<=distr_axis<=Nn
	
	w = workers()    

	wmax = min(length(w), array_size[Na+distr_axis])

	return (array_size, w[1:wmax], setindex!(d, wmax, Na+distr_axis))

end 



function _init_storage(
											 s::Tuple{Vararg{Int}},
											 ns::Tuple{Vararg{Int}},
											T::DataType;
											parallel::Bool=false,
											shared::Bool=false,
											kwargs...
											)::AbstractArray{<:Number}

	array_size, array_workers, array_distrib = _prep_init_array(s,ns; 
																						parallel=parallel, 
																						kwargs...)

	if parallel 
		
		if shared 
			
			return SharedArray{T}(array_size; pids=union(myid(), array_workers)) 

		else  # distributed array 

			return dzeros(T, array_size, array_workers, array_distrib)

		end 

	else 

		shared && @warn "SharedArray not initialized when parallel=false"

		return zeros(T, array_size)

	end 

end 

function _inds_distrib_array(array_size::NTuple{N,Int},
														array_workers::AbstractVector{Int},
														array_distrib::AbstractVector{Int};
														kwargs...
														)::Dict{Int, NTuple{N,UnitRange{Int}}
																		 } where N

	
	parallel = check_parallel_possible(;kwargs...)

	I = findall(>(1), array_distrib)

	@assert length(I)==parallel 

	t0 = Tuple(1:a for a=array_size)

	parallel || return Dict(only(array_workers)=>t0)

	i = only(I)

	distr_i = Utils.EqualDistributeBallsToBoxes_cumulRanges(
																array_size[i], array_distrib[i])
	
	return Dict(w=>Tuple(i==d ? J : a for (d,a)=enumerate(t0)) for (w,J)=zip(array_workers,distr_i))

end 


function _inds_distrib_array(
														a...; k...
														)::Dict{Int,Tuple{Vararg{UnitRange{Int}}}}

	_inds_distrib_array(_prep_init_array(a...; k...)...; k...)

end 


function inds_distrib_array(args...; 
														custom_ak::Function=init_storage_ak,
														kwargs...
														)::Dict{Int,Tuple{Vararg{UnitRange{Int}}}}

	wrapper_ak(custom_ak, _inds_distrib_array, args...; kwargs...) 

end 

#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


"""
Prepare args and kwargs for the internal method _init_storage 
"""
function init_storage_ak(T::DataType, ns::Int...; kwargs...)

	init_storage_ak(T, (), ns...)

end 

#confused with init_storage_ak(item1::Int, n::Int, kwargs...)
#function init_storage_ak(ns::Int...; kwargs...)
#	
#	init_storage_ak((), ns; kwargs...)
#
#end 
function init_storage_ak(n::Int; kwargs...)

	init_storage_ak((), (n,); kwargs...)

end 


function init_storage_ak(
												 s::Tuple{Vararg{Int}}, 
												 n::Vararg{Int}; 
												 kwargs...)

	init_storage_ak(s, n; kwargs...)

end 


function init_storage_ak(T::DataType, 
												 s::Tuple{Vararg{Int}}, 
												 n::Vararg{Int}; 
												 kwargs...)

	init_storage_ak(T, s, n; kwargs...)

end 


function init_storage_ak(T::DataType, 
												 s::Tuple{Vararg{Int}}, 
												 n::Tuple{Vararg{Int}}; 
												 kwargs...)

	(s, n, T), kwargs

end 



function init_storage_ak(
												 s::Tuple{Vararg{Int}}, 
												 n::Tuple{Vararg{Int}},
												 T::DataType=Float64;
												 kwargs...)

	(s, n, T), kwargs

end 


function init_storage_ak(item1::AbstractArray{T}, n::Int; kwargs...) where T

	init_storage_ak(item1, (n,); kwargs...) # differs from WLO 

end 


function init_storage_ak(item1::AbstractArray{T}, 
												 ns::Tuple{Vararg{Int}}; 
												 kwargs...) where T<:Number

	init_storage_ak(size(item1), ns, T; kwargs...)

end 

function init_storage_ak(item1::Number, args; kwargs...) 

	init_storage_ak([item1], args...; kwargs...)

end 



"""
Get args and kwargs needed to initialize similar storage.

"""
function similar_storage_ak(size_array_mesh::NTuple{N,Int},
												 T::DataType...;
												 kwargs...
												) where N

	@assert N>=2 
	@warn "only hotfix"

	(size_array_mesh[1:2], 
	 Tuple(nr_kPoints_from_mesh(size_array_mesh,d,N-2) for d=1:N-2),
	 T...
	 ), kwargs


end 





#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#




"""
Prepare args and kwargs using `init_storage_ak` and initiate zeros 
"""
function init_storage(args...; kwargs...)

	wrapper_ak(init_storage_ak, _init_storage, args...; kwargs...)

end  

function init_storage(item1::Dict, args...; kwargs...)::Dict

	Dict(k=>init_storage(v, args...; kwargs...) for (k,v)=item1)

end 







"""
Obtain both `init_storage` and `inds_distrib_array` in one go 

See also [`init_storage`](@ref), [`inds_distrib_array`](@ref).
"""
function init_storage_ida(args...; 
													custom_ak::Function=init_storage_ak,
													kwargs...)

	wrapper_ak(custom_ak, (_init_storage, _inds_distrib_array), args...; kwargs...)

end   



function init_storage_ida(item1::AbstractVector{<:AbstractArray}, 
													args...; kwargs...)

	AKs = [init_storage_ak(it, args...; kwargs...) for it=item1]

	return ([_init_storage(a...; k...) for (a,k)=AKs], 
					[_inds_distrib_array(a...; k...) for (a,k)=AKs]
					) 

end   



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function select_mesh_point(data::AbstractArray{T,Ns},
													 i::NTuple{Nm,Int} 
													 )::AbstractArray{T} where {T<:Number,Ns,Nm}

	view(data, fill(:,Ns-Nm)..., i...)

end 


#select_mesh_point(data, i::NTuple{N<Nm}
#									Nm, dir::NTuple{N} 


#function get_mesh_dim(Ns::Int, Nm::Int, dir::Int)::Int 
#
#	@assert Ns>=Nm && 1<=dir<=Nm 
#
#	Ns-Nm+dir 
#
#end 

#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#




"""
Internal method for in-place mesh storing. 

Adapted from WLO 
""" 
function _store_on_mesh!(f!::Function, 
												 dest::T,
												 inds::Tuple{Vararg{AbstractVector{Int}}},
												data...; 
												parallel::Bool=false,
												array_distrib=nothing,
												kwargs...
												)::T where T<:Union{<:SubOrArray{<:Number},
																		<:SubOrSArray{<:Number},
																		<:AbstractVector{<:AbstractArray},
#																		<:AbstractDict, ? 
																		}

	for i in Base.product(inds...)

#		store_on_mesh_one!(get_one!, dest, i, data...; kwargs...)
		f!(dest, i, data...; kwargs...)

	end 

	return dest 

end    


#function store_on_mesh_one!(
#							get_one_item!::Function, 
#							dest::T,
#							i::Tuple{Vararg{Int}},
#							data...; kwargs...
#												)::T where T<:Union{
#													<:SubOrArray{<:Number},
#													<:SubOrSArray{<:Number},
#													<:AbstractVector{<:AbstractArray},
#													<:AbstractDict, # ?
#													}
#
#	get_one!(select_mesh_point(dest, i), i, data...; kwargs...)
#
#	return dest 
#
#end  
#


#
#function store_on_mesh!(f!::Function, 
#												dest,
#												ns::Union{Int,Tuple{Vararg{Int}}}, 
#												args...; kwargs...)
#
#	store_on_mesh!(f!, dest, mesh_axes(ns; kwargs...), args...; kwargs...)
#
#end 
#
#
#function store_on_mesh!(f!::Function, dest::T,
#												 inds::Tuple{Vararg{AbstractVector{Int}}},
#												data...; kwargs...
#												)::T where T<:Union{<:SubOrArray{<:Number},
##																		<:SubOrSArray{<:Number}, # different
#																		<:AbstractVector{<:AbstractArray},
##																		<:AbstractDict,
#																		}
#
#	_store_on_mesh(f!, dest, inds, data...; kwargs...)
#
#end 

	#for d in data 

	#	@assert !isa(d,SubOrDArray)

	#	isa(d,SubOrSArray) && continue 
	#
	#	@assert !isa(d, SubOrDArray)

	#	isa(d,SubOrArray) || continue 
	#	
	#	length(d) < max(100, sum(length, dest)/50) && continue 

	#	@warn "Large non-shared array passed? $(size(d))"

	#	@show typeof(d)
	#end 


#function store_on_mesh!(f!::Function, dest::T, 
#												 inds::Tuple{Vararg{AbstractVector{Int}}},
#												data...; 
#												array_distrib::Dict{Int,Tuple{Vararg{UnitRange{Int}}}},
#												kwargs...
#												)::T where T<:Union{<:AbstractVector{<:SubOrSArray},
#																		 SubOrSArray{<:Number},
#																		 }
#
#@warn "Not implemented" 
#
##	map(procs_(array_distrib)) do p 
##
##		locinds = array_distrib[p]
##
##		@spawnat p begin   
##
##			_store_on_mesh!(f!, 
##											localpart_(dest, locinds), # NOT DEFINED 
##											map(findall_indexin, get_big_inds(locinds), inds), 
###											pass locinds  ? 
##										data...;
##										kwargs...)
##
##		end  
##
##	end .|> fetch 
#
#	return dest 
#
#end   




















































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































############################################################################ 
end # module MeshInPlace
