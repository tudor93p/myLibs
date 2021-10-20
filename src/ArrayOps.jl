module ArrayOps
#############################################################################

import ..LA, ..SpA 

import ..Utils




#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#




function init_zeros(A::Union{AbstractArray,Tuple},
										d_to_n::Vararg{Pair{Int,Int}},
										)::Array

	init_zeros(rand(A), d_to_n...)

end  

function init_zeros(x::T,
										d_to_n::Vararg{Pair{Int,Int},N},
										)::Array{T,N} where {T<:Number,N}

	init_zeros(T, d_to_n...)

end 


function init_zeros(T::DataType,
										d_to_n::Vararg{Pair{Int,Int},N},
										)::Array{T,N} where N

	shape = zeros(Int, N) 

	for (d,n) in d_to_n 

		shape[d] = n 

	end 

	return zeros(T, shape...)

end 



function init_zeros(d_to_n::Vararg{Pair{Int,Int},N}
									 )::Array{Float64,N} where N 
	
	init_zeros(Float64, d_to_n...)

end



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function findLostDims(A::AbstractArray, B::Number, args...)::Vector{Int}

	1:ndims(A)

end


function findLostDims(A::AbstractArray, B::AbstractArray, possible_dims=[])::Vector{Int}


	@assert ndims(A)>=ndims(B) "B cannot be subarray of A"


	#root((sA,sB,p))::Dict{Int,Int} = Dict{Int,Int}()


	function possible_extensions((sA,sB,p), candidate::Dict)::Vector{Dict}
		
		n = length(candidate)+1

		n>length(sB) && return []

		length(sA)==length(sB) && return [Dict(i=>i for i in 1:length(sA))]

		start = max(maximum(values(candidate),init=0)+1,n)

		stop = n + length(sA)-length(sB)

		return [merge(candidate, Dict(n=>i)) for i in start:stop]

	end 
		


	function promising_candidate((sA,sB,p), candidate::Dict)::Bool

		if !isempty(candidate) 

			n = length(candidate)

			if length(sA)-length(sB)>0 && !isempty(p) 

				count(in(values(candidate)),p)>length(sA)-length(sB) && return false

			end 

			in(sB[n], [1, sA[candidate[n]]]) || return false 

			n<2 || candidate[n]>candidate[n-1] || return false 

		end

		return true 
			
	end 



	function accept_sol((sA,sB,p), candidate::Dict)::Bool

		length(candidate)==length(sB) && issubset(values(candidate),1:length(sA))

	end 

	solutions = Utils.Backtracking(
							 (size(A),size(B),vcat(possible_dims...)), 
							 	possible_extensions,
							 	promising_candidate,
							 	accept_sol,
								x->Dict(),
							 	)

	out = [setdiff(1:ndims(A), values(s)) for s in solutions]

	length(out)==1 && return out[1]

	length(out)==0 && error("No solution was found. Check the sizes: ",size(A)," ",size(B))

	error("\n",out,"\n'A' has the same size along multiple dimensions and the lost dimensions cannot be determined uniquely" )

end 



function restoreLostDims(f::Function, A::AbstractArray, 
												 args...)::AbstractArray

	restoreLostDims(A, f(A), args...)

end

function restoreLostDims(A::AbstractArray, B::AbstractArray,
												 args...)::AbstractArray

	lost_dims = findLostDims(A, B, args...)

	isempty(lost_dims) && return B

	new_size = (in(i,lost_dims) ? 1 : size(A,i) for i=1:ndims(A))

	return reshape(B, new_size...)

end 


#lost_dims = findLostDims(A, B, dims) 




function newPosIndex_afterLostDims(A::AbstractArray, B::AbstractArray, 
																	 args...)::Function
	
	lost_dims = findLostDims(A, B, args...)

	return n -> n - count(<(n), lost_dims)

end 

function newPosIndex_afterLostDims(n::Int, args...)::Int

	newPosIndex_afterLostDims(args...)(n)

end 







#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function ApplyF_IndsSets(f::Function, 
												 A::AbstractArray{<:Number},
												 slice_dim::Int, 
												 inds::AbstractVector{<:AbstractVector{Int}};
												 kwargs...
												 )::AbstractVector{<:AbstractArray{<:Number}}

	ApplyF_IndsSets(f, A, slice_dim, inds, setdiff(1:ndims(A),slice_dim);
									kwargs...)

end 




function ApplyF_IndsSets(f::Function, 
												 A::AbstractArray, 
												 slice_dim::Int, 
												 inds::AbstractVector{<:AbstractVector{<:Int}},
												 f_dim::Union{Int,AbstractVector{Int}};
												 keepdims=false,
												 )::AbstractVector{AbstractArray}

	isempty(intersect(slice_dim, f_dim)) || error("Axes must be disjunct")

	a = selectdim(A, slice_dim, vcat(inds...))

	if keepdims

		return Slice_LikeArraysInList(inds, 
																	restoreLostDims(f, a, f_dim), 
																	slice_dim)
	end 

	b = f(a)

	return Slice_LikeArraysInList(
								inds, b, 
								newPosIndex_afterLostDims(slice_dim, a, b, f_dim)
												 )

end



function ApplyF_IndsSets_(
													f::Function, 
												  A::AbstractArray, 
												  slice_dim::Int, 
												  inds::AbstractVector{<:AbstractVector{<:Int}},
												  f_dim::Union{Int,AbstractVector{Int}};
												  keepdims=false,
												  )::AbstractVector{AbstractArray}

	map(inds) do i

		a = selectdim(A, slice_dim, i) 
	
		return keepdims ? restoreLostDims(f, a, f_dim) : f(a)

	end 

end 


function ApplyF_IndsSets_Check(args...; kwargs...)


	result1 = ApplyF_IndsSets(args...; kwargs...)

	@time result1= ApplyF_IndsSets(args...; kwargs...)
	
	result2 = ApplyF_IndsSets_(args...; kwargs...)

	@time result2 = ApplyF_IndsSets_(args...; kwargs...)


	for p in zip(result1,result2)

		if isapprox(p...) 
			
			println("test passed")

		else 

			error("test not passed")

		end 

	end

	return result1


end 






#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#



function mapslices_dropLostDims(f::Function, A::AbstractArray, dims)

	dropLostDims(A, mapslices(f, A, dims=dims), dims)

end 

function dropLostDims(f::Function, A::AbstractArray, dims)

	dropLostDims(A, f(A), dims)

end 

function dropLostDims(A::AbstractArray, B::AbstractArray, dims)

	not_lost_dims = setdiff(dims, findLostDims(A, B, dims))

	isempty(not_lost_dims) && return B

	return dropdims(B, dims=Tuple(not_lost_dims))

end 


function softDropdims(A::AbstractArray, dims)

	D = intersect(findall(size(A).==1),dims)

	isempty(D) && return A

	return dropdims(A, dims=Tuple(D))

end 

#function softDropdims(f::Function, A::AbstractArray, dims)
#
#	B = f(A) 
#
#	lost_dims = findLostDims(A, B, dims)
#	
#end 



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

multiply_elementwise(A::Number, B::AbstractArray)::Array = A*B
multiply_elementwise(A::AbstractArray, B::Number)::Array = A*B
multiply_elementwise(A::Number, B::Number)::Number = A*B

function multiply_elementwise(A::AbstractArray, B::AbstractArray, dim=nothing)::Array

	size(A)==size(B) && return A.*B

	sA,sB = filter.(!isequal(1), size.((A,B)))

	inds = filter(!isnothing,  indexin(vcat(sA...), vcat(sB...)))

	if length(inds)!=minimum(length,(sA,sB)) || any(diff(inds).<=0)

		error("Shape mismatch")
	
	end 

	ndims(A)==ndims(B) && return A.*B


	if issubset(sA, sB) # A is smaller than B => expand A

		new_size_A,inds = vcat(size(B)...), vcat(1:ndims(B)...)

		for i in inds 

			!isnothing(dim) && i==dim && continue

			nr = count(new_size_A[i].==sA)

			nr==0 && setindex!(new_size_A, 1, i)

			nr>1 && error("Ambiguous result. Too many occurences found")

		end

		return reshape(A, new_size_A...).*B
		

	elseif issubset(sB, sA)
		
		return multiply_elementwise(B, A)

	end 

	error()

end 








#===========================================================================#
#
#	Dense matrix from lists of indices and values
#				(each row in inds represents a catesian index)
#
#---------------------------------------------------------------------------#


function Array_from_ListIndsVals(inds::AbstractMatrix{Int}, 
																 vals::AbstractVector{T}, 
																 shape::AbstractVector{Int}
																 )::Array{T} where T<:Number

	v0 = first(vals)

	i0 = fill(:, ndims(v0))

	A = zeros(typeof(first(v0)), size(v0)..., shape...)

	for (I,V) in zip(eachrow(inds), vals)
		
		A[i0..., I...] = V 

	end 

	return A

end



function Array_from_ListIndsVals(inds::Tuple, args...)::Array

	Array_from_ListIndsVals(collect(inds), args...)

end

function Array_from_ListIndsVals(inds::AbstractVector{<:Number}, 
																 args...)::Array

	Array_from_ListIndsVals(Utils.VecAsMat(inds, 2), args...) # column vector 

end 


function Array_from_ListIndsVals(inds::AbstractVector{<:AbstractVector{<:Number}}, args...)::Array

	Array_from_ListIndsVals(vcat((Utils.VecAsMat(i, 1) for i in inds)...),
													args...)

end 



function Array_from_ListIndsVals(inds::AbstractMatrix, vals)::Array

	Array_from_ListIndsVals(inds, vals, maximum.(eachcol(inds)))

end 


#function Array_from_ListIndsVals(inds::AbstractMatrix, vals::Utils.List, args...)::Array  
#
#	Array_from_ListIndsVals(inds, collect(vals), args...)
#
#end 



#===========================================================================#
#
#	Cumulative sum for a list of lists
#
#---------------------------------------------------------------------------#

function recursive_end(array::T) where T


	T <: Number && return array

	isempty(array) && return 0

	T <: AbstractArray{<:Number} && return array[end]

	i = findlast(!isempty,array)

	isnothing(i) && return 0

	return recursive_end(array[i])


end



function recursive_cumsum(array::T,start=0) where T

	T <: Number && return array + start

	nonempty_arrayinds = findall(!isempty,array)


	isempty(nonempty_arrayinds) && return array


	T <: AbstractArray{<:Number} && return cumsum(array) .+ start


	out = copy(array)

	for (order_nonempty,array_index) in enumerate(nonempty_arrayinds)

		out[array_index] = recursive_cumsum(array[array_index],

																				if order_nonempty==1 start else

									recursive_end(out[nonempty_arrayinds[order_nonempty-1]]) 

																				end)
	end

	return	out

end







#===========================================================================#
#
#	Block diagonal matrix 
#
#---------------------------------------------------------------------------#

function BlkDiag(arg::T)::Matrix where T
	
	Utils.isList(T,AbstractMatrix) && return cat(arg...,dims=(1,2))

	T<:Base.Generator &&  return cat(arg...,dims=(1,2))

	T<:AbstractMatrix{<:Number} && return arg

	T<:AbstractVector{<:Number} && return LA.diagm(0=>arg)

	T<:Number && return hcat(arg)

	Utils.isList(T,Any) && return BlkDiag(map(BlkDiag,arg))

	error("Input not understood. Type: ",T)

end

BlkDiag(args...) = BlkDiag(args)







#===========================================================================#
#
# Unit matrix (built-in expression not intuitive)
#
#---------------------------------------------------------------------------#

function UnitMatrix(d::Int64, dtype::DataType=Float64)::Matrix

  Matrix{dtype}(LA.I,d,d)

end

function UnitMatrix(A::AbstractMatrix)::Matrix

	one(A)

end

function UnitMatrix(n::Nothing=nothing)::Nothing 

	nothing

end








#===========================================================================#
#
# Make a big sparse array with the elements of given array at certain indices
#
#---------------------------------------------------------------------------#

function Bigger_Matrix(new_Is,new_Js,size_Is,size_Js,matrix)

	# By converting matrix to sparse
    I,J,V = SpA.findnz(SpA.sparse(matrix))

	# Withouth converting to sparse -- much slower!
#    IJ = findall(!iszero,matrix)

#    (I,J),V = vcat.(Tuple.(IJ)...), getindex.([matrix],IJ)

    return SpA.sparse(new_Is[I],new_Js[J],V,size_Is,size_Js)

end



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#



function Slice_IndsSets(inds::AbstractVector{<:AbstractVector{Int}}, 
												X::AbstractArray, dim::Int=1)::Vector

	Slice_LikeArraysInList(inds, selectdim(X, dim, vcat(inds...)), dim)

end






function Slice_LikeArraysInList(arrays, X::AbstractArray, dim::Int=1
																)::Vector{<:Array}

	[collect(selectdim(X, dim, m)) for m in Utils.sepLengths_cumulRanges(arrays)]

end







#===========================================================================#
#
# Transforms a n-dimensional array intro a list of (n-1)-dimensional arrays
#
#---------------------------------------------------------------------------#


function ArrayToList(a)

  return [a[i,axes(a)[2:ndims(a)]...] for i in axes(a,1)]

end





#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

#===========================================================================#
#
# Normalize an array (full, on columns, on rows)
#
#---------------------------------------------------------------------------#

function Normalize(A::AbstractArray{T, n} where T, p::Real=2; 
									 tol::Float64=1e-12, kwargs...)::Array{<:Number, n} where n

	if n>1 && haskey(kwargs, :dim) 
	
		dim = kwargs[:dim]

		@assert isa(dim,Int) && 1<=dim<=n 

		N = LA.norm.(eachslice(A, dims=dim), p)

		N[N.<tol] .= tol

		return A./reshape(N, (d==dim ? Colon() : 1 for d=1:n)...)

	end 

	return A/LA.norm(A, p)

end


function Normalize_Columns(A::AbstractMatrix{<:Number}, p...)::Matrix 
	
	Normalize(A, p...; dim=2)

end 

function Normalize_Rows(A::AbstractMatrix{<:Number}, p...)::Matrix 
	
	Normalize(A, p...; dim=1)

end 


























































































































































































































































































































































#############################################################################
end 

