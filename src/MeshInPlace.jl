module MeshInPlace
############################################################################

import Distributed: @spawnat, myid, workers, nworkers, procs
#import DistributedArrays: DArray, SubOrDArray, localindices, localpart, dzeros 

using SharedArrays:SharedArray


import ..Utils 



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


"""
Size of a storage array. 

# Arguments
- `s::NTuple{Na,Int}`: size of the array `A` at each mesh point, `ndims(A)==Na` 
- `ns::NTuple{Nm,Int}`: mesh size, dimensionality `Nm`

# Keyword arguments 
- `endpoint::Bool=true`: whether the last and first mesh points differ
"""
function _storage_size(
											s::NTuple{Na,Int},
											ns::NTuple{Nm,Int},
											args...;
											endpoint::Bool=true,
											kwargs...
											)::NTuple{Na+Nm,Int} where {Na,Nm}

	Utils.tuplejoin(s,endpoint ? ns : ns .- 1)

end 





#===========================================================================#
#
# _prep_init_array and _init_storage take identical same args and kwargs
# internal methods 
#
#---------------------------------------------------------------------------#


function _prep_init_array(
											s::NTuple{Na,Int},
											ns::NTuple{Nn,Int},
											args...;
											distr_axis::Int=Nn,
											kwargs...
											)::Tuple{NTuple{Na+Nn,Int}, Vector{Int}, Vector{Int}
															 } where {Na,Nn}

	array_size = _storage_size(s, ns; kwargs...)

	d = ones(Int, Na+Nn)


	parallel = check_parallel_possible(;kwargs...)


	parallel || return (array_size, [myid()], d) # just the master proc

	@assert 1<=distr_axis<=Nn
	
	w = workers()    

	wmax = min(length(w), array_size[Na+distr_axis])

	return (array_size, w[1:wmax], setindex!(d, wmax, Na+distr_axis))

end 



function _init_storage(
											s::NTuple{Ns,Int},
											ns::NTuple{Nn,Int},
											T::DataType;
											parallel::Bool=false,
											shared::Bool=false,
											kwargs...
											)::AbstractArray{<:Number,Ns+Nn} where Ns where Nn 

	array_size, array_workers, array_distrib = _prep_init_array(s,ns; 
																						parallel=parallel, 
																						kwargs...)


	parallel || return zeros(T, array_size)

	shared || return dzeros(T, array_size, array_workers, array_distrib)

	return SharedArray{T}(array_size; pids=union(myid(), array_workers))

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
												 T::DataType;
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




























































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































############################################################################ 
end # module MeshInPlace
