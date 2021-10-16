
const TOLERANCE = 1e-5 
approx(args...) = isapprox(args...; atol=TOLERANCE)



const VECTOR_STORE_DIM = 2 
const VECTOR_AXIS = [2,1][VECTOR_STORE_DIM] 


Unique(args...;kwargs...) = Utils.Unique(args...; kwargs..., tol=TOLERANCE, dim=VECTOR_STORE_DIM)
Unique!(args...;kwargs...) = Utils.Unique!(args...; kwargs..., tol=TOLERANCE, dim=VECTOR_STORE_DIM)


EnumUnique(args...;kwargs...) = Utils.EnumUnique(args...; kwargs..., tol=TOLERANCE,dim=VECTOR_STORE_DIM)

@assert Set([VECTOR_STORE_DIM,VECTOR_AXIS])==Set(1:2)

Cat = [vcat,hcat][VECTOR_STORE_DIM]

eachvec = [eachrow,eachcol][VECTOR_STORE_DIM]


function VecsOnDim(vecs::T; dim::Int)::T where T<:AbstractMatrix{<:Number}

	dim==VECTOR_STORE_DIM && return vecs 

	return transpose(vecs)

end 





sortvecs(A::AbstractMatrix{<:Number}; kwargs...) = sortslices(A; dims=VECTOR_STORE_DIM, kwargs...)

VecsToMat(args...; kwargs...) = Utils.VecsToMat(args...; dim=VECTOR_STORE_DIM, kwargs...)

vectors_of_integers(args...; kwargs...) = Utils.vectors_of_integers(args...; dim=VECTOR_STORE_DIM, kwargs...) 



CombsOfVecs10(args...; kwargs...) = Utils.CombsOfVecs10(args...; dim=VECTOR_STORE_DIM, kwargs...)

PointInPolygon(args...; kwargs...) = Geometry.PointInPolygon_wn(args...; dim=VECTOR_STORE_DIM, kwargs...)

#===========================================================================#
#
# Number of vectors and the indices of the vectors in a matrix
#
#---------------------------------------------------------------------------#


NrVecs(A::AbstractMatrix{<:Number}, x::Nothing=nothing)::Int = size(A, VECTOR_STORE_DIM)   

IndsVecs(A::AbstractMatrix{<:Number})::Vector{Int} = axes(A, VECTOR_STORE_DIM)


function IndsVecs(I::AbstractVector{Int}, d::AbstractArray{Int, N}
									)::Array{Int, N} where N 

	@assert issubset(d, I) "Index $d not possible"
	
	return d

end 

function IndsVecs(L::Int,
									d::AbstractVector{Bool})::Vector{Bool}

	@assert L==length(d) "Index $d not possible"
	
	return d

end 


function IndsVecs(I::AbstractVector{Int}, d::Int)::Vector{Int}

	IndsVecs(I, [d])

end 


function IndsVecs(A::AbstractMatrix{<:Number}, 
									d::Union{Int, AbstractArray{Int}})::Array{Int}
	
	IndsVecs(IndsVecs(A), d)

end  

function IndsVecs(A::AbstractMatrix{<:Number}, 
									d::AbstractVector{Bool})::Vector{Bool}
	
	IndsVecs(NrVecs(A), d)

end  



function NrVecs(A::AbstractMatrix{<:Number}, 
								d::Union{Int, AbstractVector{Int}})::Int

	length(IndsVecs(A, d))
	
end 



#===========================================================================#
#
# Get the actual vectors
#
#---------------------------------------------------------------------------#


function Vecs(A::AbstractMatrix{T})::Matrix{T} where T<:Number 
	
	A

end 

function Vecs(A::AbstractMatrix{T}, 
							d::AbstractVector{Bool})::Matrix{T} where T<:Number 

	selectdim(A, VECTOR_STORE_DIM, IndsVecs(A, d))

end 

function Vecs(A::AbstractMatrix{T}, d::AbstractArray{Int, N}
								 )::Array{T, N+1} where T<:Number where N

	selectdim(A, VECTOR_STORE_DIM, IndsVecs(A, d))

end 


function Vecs(A::AbstractMatrix{T}, d::Int)::Matrix{T} where T<:Number 

	Vecs(A, [d])

end 



function SetVecs!(A::AbstractMatrix{T}, inds_vecs, X, Xinds=Colon()
								 )::Matrix{T} where T<:Number
	
	setindex!(selectdim(A, VECTOR_STORE_DIM, inds_vecs), X, Xinds)
	
end 


#===========================================================================#
#
# Length and axis of individual vectors 
#
#---------------------------------------------------------------------------#


VecLen(A::AbstractMatrix{<:Number})::Int = size(A, VECTOR_AXIS)
VecLen(v::AbstractVector{<:Number})::Int = length(v)

VecAx(A::AbstractMatrix{<:Number})::Vector{Int} = axes(A,VECTOR_AXIS)
VecAx(A::AbstractVector{<:Number})::Vector{Int} = 1:length(v)





#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

Shape()::Vector{Int} = zeros(Int,2)

function Shape(vec_len::Int)::Vector{Int}

	setindex!(Shape(), vec_len, VECTOR_AXIS) 

end 

function Shape(vec_len::Int, nr_vecs::Int)::Vector{Int}

	setindex!(Shape(vec_len), nr_vecs, VECTOR_STORE_DIM)

end 



function ZeroVecs(args...)::Matrix{Float64}

	zeros(Float64, Shape(args...)...)

end 



function EmptyVec(vec_len::Int64)::Matrix{Float64}

	ZeroVecs(vec_len) 

end 

EmptyVec(n::Nothing=nothing)::Matrix{Float64} = EmptyVec(0)
 

function EmptyVec(M::AbstractMatrix{Float64})::Matrix{Float64}
	
	EmptyVec(VecLen(M))

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

function to_myMatrix(x::Nothing, D::Union{Int64,Nothing}
										)::Matrix{Float64}

	EmptyPos(D)

end 


function to_myMatrix(x::Real, D::Union{Int64,Nothing}
										)::Matrix{Float64} 

	fill(x, Shape(Utils.Assign_Value(D, 1), 1)...)

end 



function to_myMatrix(x::AbstractMatrix{<:Number}, D::Union{Int64,Nothing}=nothing
										)::Matrix{Float64}

	@assert isnothing(D) || D==VecLen(x) "Wrong array dimension" 

	return x
	
end 


function to_myMatrix(x::Utils.List, D::Union{Int64,Nothing}=nothing
										)::Matrix{Float64}

	isempty(x) && return EmptyPos(D) 


	if Utils.isList(x, Real) && (isnothing(D) || length(x)==D)

		return Utils.VecAsMat(collect(x), VECTOR_STORE_DIM)


	elseif all(Utils.isList(;T=Real), x) && (
					 isnothing(D) || all(D.==length.(x)))

		out = ZeroVecs(length(x[1]), length(x))
	
		for (j,xj) in enumerate(x) # each column j is a vector

			SetVecs!(out, j, collect(xj))

		end 

		return out 

	end 

	println()

	println.(typeof.(x))

	println()
	println.(x)
	println()
	error("Wrong dimensions or types")

end 



function to_myMatrix(x, M::AbstractMatrix{<:Number})::Matrix{Float64}

	to_myMatrix(x, VecLen(M))
	
end 


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function CombsOfVecs(A::AbstractMatrix{Float64}, coeff)::Matrix{Float64}

	Utils.CombsOfVecs(A, to_myMatrix(coeff, NrVecs(A)); dim=2)
	
end 


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#



function to_myODict(x::Nothing=nothing, D::Union{Int64,Nothing}=nothing
									 )::OrderedDict{Any, Matrix{Float64}}

	OrderedDict()

end


function to_myODict(x::AbstractDict, D::Union{Int64,Nothing}=nothing
									 )::OrderedDict{Any, Matrix{Float64}}

	if x isa OrderedDict 

		valtype(x)<:AbstractMatrix{<:Float64} && return x

		return OrderedDict(k=>to_myMatrix(v, D) for (k,v) in x)

	end 

	return OrderedDict(map(sort(collect(keys(x)), by=string)) do k

									 k => to_myMatrix(x[k], D)

								 end)

end 






function to_myODict((labels,atoms)::Tuple{<:AbstractVector{<:AbstractString},
																					<:Any},
										D::Union{Int64,Nothing}=nothing 
										)::OrderedDict{Any, Matrix{Float64}}

	to_myODict((labels,to_myMatrix(atoms, D)),D)

end 

function to_myODict((labels,atoms)::Tuple{<:AbstractVector{<:AbstractString},
																					<:AbstractMatrix{<:Number}},
										D::Union{Int64,Nothing}=nothing 
										)::OrderedDict{Any, Matrix{Float64}}

	l,a = length(labels), NrVecs(atoms)

	if a>l && a%l==0

		return to_myODict((repeat(labels,outer=div(a,l)),atoms), D)

	elseif a!=l 

		error("Sizes don't match")

	end 


	return OrderedDict(L=>Vecs(atoms, I) for (L,I) in EnumUnique(labels))



end 






function to_myODict(x::Pair, D::Union{Int64,Nothing}=nothing
									 )::OrderedDict{Any, Matrix{Float64}}

	OrderedDict(x.first => to_myMatrix(x.second, D))

end


function to_myODict(x::AbstractMatrix{<:Number}, D::Union{Int64,Nothing}=nothing
									 )::OrderedDict{Any, Matrix{Float64}}

	OrderedDict("A"=>to_myMatrix(x, D))

end



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function to_myODict(x::Utils.List, D::Union{Int64,Nothing}=nothing
									 )::OrderedDict{Any, Matrix{Float64}}

	isempty(x) && return OrderedDict()


	types = map(x) do xi 

		if xi isa Union{<:Pair,
										<:AbstractDict, 
										<:Tuple{<:AbstractVector{<:AbstractString},	<:Any}}
			return 1 

		elseif xi isa AbstractMatrix{<:Number}

			return 2

		else 

			return 0

		end

	end 


	all(types.==0) && return to_myODict(to_myMatrix(x,D), D)
#it's just one matrix given in a strange way



	@assert all(types.>0) "Some input not understood"


	out1 = merge(Cat, to_myODict(), to_myODict.(x[types.==1], D)...) 


	function good_keys(used_keys, n::Int, sol=[])
		
		length(sol)==n && return sol 

		start = isempty(sol) ? 'A' : sol[end][1]+1

		for K in zip(lowercase(start):'z', uppercase(start):'Z')

			any(s(k) in used_keys 
						for k in K for s in (string, Symbol)) && continue

			return good_keys(used_keys, n, [sol;string(K[2])])

		end 

		error("Not enough trials")

	end 

	unused_keys = good_keys(keys(out1), count(types.==2))


	return merge!(Cat, out1, 
								(OrderedDict(k=>to_myMatrix(v, D)) 
								 					for (k,v) in zip(unused_keys,x[types.==2]))...)


end





#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#






function to_myODict(x, M::AbstractMatrix{<:Number}
									 )::OrderedDict{Any, Matrix{Float64}}

	to_myODict(x, VecLen(M))

end 

 

#to_myODict(x...)::OrderedDict{Any, AbstractMatrix{Float64}} = to_myODict(x)





