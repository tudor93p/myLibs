module Lattices
#############################################################################

using OrderedCollections:OrderedDict

import ..LA

import ..Utils, ..Algebra, ..ArrayOps, ..Geometry

import Plots, ConcaveHull, QHull




#export Lattice



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

const TOLERANCE = 1e-5

approx(args...) = isapprox(args...; atol=TOLERANCE)

Unique(args...;kwargs...) = Utils.Unique(args...; kwargs..., tol=TOLERANCE)
Unique!(args...;kwargs...) = Utils.Unique!(args...; kwargs..., tol=TOLERANCE)

# Latt : python object
# Properties:
#		Latt.LattVec
# Methods which modify: 
# 	Latt.Add_atoms(sublattices=nothing)
# 	Latt.Shift_Atoms(...)
# 	Latt.Reduce_Dim(...)
# Methods which don't:
#		Latt.Distances(...)
# 	Latt.Supercell(...)
# 	Latt.PosAtoms(sublattices)
#

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

Base.copy(s::Symbol) = s
Base.copy(s::String) = s

#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


mutable struct Lattice

	LattVec::Matrix{Float64}

	LattDims::Vector{Int}

	Atoms::OrderedDict{Any, Matrix{Float64}}

	Vacancies::OrderedDict{Any, Matrix{Float64}}




	function Lattice(LV, SL=nothing, VA=nothing, LD=nothing; mode::Symbol=:cartesian)

		lv = to_myMatrix(LV) 
		
		f = Dict(:cartesian=>identity,:fractional=>transpose)[mode]

		sl, va = to_myODict(SL, f(lv)), to_myODict(VA, f(lv)) 


		latt_dim, vec_dim = LattDim(lv, LD), VecDim(lv) 




		for k in [keys(sl);keys(va)], odict in [sl,va]

			if haskey(odict, k)

				s = VecDim(odict[k])

				if mode==:fractional

					s==latt_dim || error("Dimensionality mismatch. $s coefficient(s) provided for sublattice '$k', but $latt_dim are necessary")

					odict[k] = CombsOfVecs(LattVec(lv, LD), odict[k])
				
				elseif mode==:cartesian

					vec_dim!=s && error("Dimensionality mismatch. The positions for sublattice '$k' are $s-dimensional, but the lattice vectors are $vec_dim-dimensional")

				else 

					error("Mode $mode not defined")

				end 


			end 

		end 



		return new(lv, Utils.Assign_Value(LD, 1:latt_dim), sl, va)

	end 


	function Lattice(latt::Lattice; act_on_vectors=identity,
									 								act_on_atoms = identity,
																	)

		vectors,dims = applyOnLattVec(act_on_vectors, latt)

		return Lattice(vectors, applyOnLattAtoms(act_on_atoms, latt)..., dims)

	end 
	
		
end 



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

#function Lattice!(latt::Lattice; act_on_vectors! =nothing,
#																	act_on_atoms! =nothing,	)
#
#
#	if hasmethod(act_on_vectors!, (AbstractMatrix,))
#	
#		act_on_vectors!(latt.LattVec)
#	
#	elseif hasmethod(act_on_vectors!, (Lattice,))
#	
#		act_on_vectors!(latt)
#	
#	elseif !isnothing(act_on_vectors!)
#	
#		error("Improper definiton of 'act_on_vectors!'")
#	
#	end
#
#
#	isnothing(act_on_atoms!) && return 
#
#	foreach([:Atoms, :Vacancies]) do kind
#
#		D = getproperty(latt, kind)
#
#		for args in ((D,), (kind,), (), (kind, D),)
#
#			!hasmethod(act_on_atoms!, typeof.(args)) && continue
#			
#			return act_on_atoms!(arg...)
#
#		end 
#
#
#		foreach(collect(keys(D))) do k 
#
#			for args in ((D[k],), (kind, D[k]), (D[k],k), (kind, D[k], k))
#
#				!hasmethod(act_on_atoms!, typeof.(args)) && continue 
#				
#				return act_on_atoms!(args...)
#
#			end 
#
#			error("Improper methods of 'act_on_atoms!'")
#
#		end 
#
#	end
#		
#end 
#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

Base.getindex(latt::Lattice, kind::Symbol) = getproperty(latt, kind)
Base.getindex(latt::Lattice, i::Int) = latt[[:Atoms,:Vacancies][i]]

Base.firstindex(latt::Lattice) = :Atoms
Base.lastindex(latt::Lattice) = :Vacancies

Base.length(latt::Lattice) = 2

function Base.iterate(latt::Lattice, state=(1,0))

	SV = [:Atoms, :Vacancies] 

	element, count = state 

	count >= length(SV) && return nothing

	return ((SV[element], latt[element]), (element + 1, count + 1))

end 



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#




function applyOnLattVec!(f::Function, latt::Lattice)::Lattice 
#	AbstractMatrix,AbstractVector}

	LattDims!(f, latt) 

	LattVec!(f, latt) 


	if applicable(f, latt) 

		v, d = f(latt)

		LattDims!(latt, d)

		LattVec!(latt, v) 

	end 

	return latt 

end 



function applyOnLattVec(f::Function, latt::Lattice)::Tuple{AbstractMatrix,AbstractVector}

#	applicable(f, latt) && return f(latt)

	(LattVec(f, latt), LattDims(f, latt))

end 





#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#



#			(Dict,)
#			(Symbol,)
#			()
#			(Symbol, Dict)
#
#			(Dict,k::Symbol)
#			(Symbol,k::Symbol)
#			(k::Symbol,)
#			(Symbol, Dict,k::Symbol)
#
#			(Matrix,)
#			(Symbol, Matrix)
#
#			(Matrix,k)
#			(Symbol, Matrix, k)

#	 ?? 	(D,k), (kind,k), (kind, D, k), ()


function applyOnLattAtoms(f::Function, latt::Lattice)::Vector

	map(latt) do (kind,D)

		for args in ((D,), (kind,), (kind, D),) 

			applicable(f, args...) && return f(map(copy, args)...)

		end 


		return map(sublatt_labels(D)) do k 

			for args in ((D[k],), (kind, D[k]), (D[k],k), (kind, D[k], k))

				applicable(f, args...) && return f(map(copy, args)...)

			end 

			error("No method of 'f' could be used")

		end 

	end
	
end 



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#



function applyOnLattAtoms!(f::Function, latt::Lattice)::Lattice 

	foreach(latt) do (kind,D)

		for args in ((D,), (kind,), (kind, D),) 

			applicable(f, args...) || continue 
			
			return setindex!(latt, f(map(copy, args)...), kind)

		end 

		foreach(sublatt_labels(D)) do k 

			for arg in ((D, k), (D[k], k), (D[k],),
									(kind, D, k), (kind, D[k], k), (kind, D[k]))

				applicable(f, arg...) && return setindex!(D, f(map(copy, arg)...), k)

			end 


#
			error("No method of 'f' was used")

		end 

	end
	
	return latt

end 



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function to_myMatrix(x::Nothing, D::Union{Int64,Nothing}
	)::AbstractMatrix{Float64}

	EmptyPos(D)

end 


function to_myMatrix(x::Real, D::Union{Int64,Nothing}
										)::AbstractMatrix{Float64} 

	fill(x, Utils.Assign_Value(D, 1), 1)

end 



function to_myMatrix(x::AbstractMatrix, D::Union{Int64,Nothing}=nothing
										)::AbstractMatrix{Float64}

	(isnothing(D) || size(x,1)==D) && return x

	error("Wrong array dimension")
	
end 


function to_myMatrix(x::Utils.List, D::Union{Int64,Nothing}=nothing
										)::AbstractMatrix{Float64}

	isempty(x) && return EmptyPos(D) 


	if Utils.isList(x, Real) && (isnothing(D) || length(x)==D)

			return reshape(collect(x),:,1)

	elseif all(Utils.isList(;T=Real), x) && (
					 isnothing(D) || all(D.==length.(x)))

		out = zeros(Float64, length(x[1]), length(x))
	
		for (j,xj) in enumerate(x) # each column j is a vector

			out[:,j] = collect(xj)

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


function to_myMatrix(x, arg::Union{<:Lattice,<:AbstractMatrix}
										)::AbstractMatrix{Float64}

	to_myMatrix(x, VecDim(arg))
	
end 


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#



function to_myODict(x::Nothing=nothing, D::Union{Int64,Nothing}=nothing
									 )::OrderedDict{Any, AbstractMatrix{Float64}}

	OrderedDict()

end


function to_myODict(x::AbstractDict, D::Union{Int64,Nothing}=nothing
									 )::OrderedDict{Any, AbstractMatrix{Float64}}

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
										)::OrderedDict{Any, AbstractMatrix{Float64}}

	to_myODict((labels,to_myMatrix(atoms, D)),D)

end 

function to_myODict((labels,atoms)::Tuple{<:AbstractVector{<:AbstractString},
																					<:AbstractMatrix},
										D::Union{Int64,Nothing}=nothing 
										)::OrderedDict{Any, AbstractMatrix{Float64}}

	l,a = length(labels), size(atoms,2)

	if a>l && a%l==0

		return to_myODict((repeat(labels,outer=div(a,l)),atoms), D)

	elseif a!=l 

		error("Sizes don't match")

	end 

	return OrderedDict(L=>atoms[:,I] 
										 for (L,I) in zip(Unique(labels, inds=:all)...))

end 






function to_myODict(x::Pair, D::Union{Int64,Nothing}=nothing
									 )::OrderedDict{Any, AbstractMatrix{Float64}}

	OrderedDict(x.first => to_myMatrix(x.second, D))

end


function to_myODict(x::AbstractMatrix, D::Union{Int64,Nothing}=nothing
									 )::OrderedDict{Any, AbstractMatrix{Float64}}

	OrderedDict("A"=>to_myMatrix(x, D))

end



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function to_myODict(x::Utils.List, D::Union{Int64,Nothing}=nothing
									 )::OrderedDict{Any, AbstractMatrix{Float64}}

	isempty(x) && return OrderedDict()


	types = map(x) do xi 

		if xi isa Union{<:Pair,
										<:AbstractDict, 
										<:Tuple{<:AbstractVector{<:AbstractString},	<:Any}}
			return 1 

		elseif xi isa AbstractMatrix

			return 2

		else 

			return 0

		end

	end 


	all(types.==0) && return to_myODict(to_myMatrix(x,D), D)
#it's just one matrix given in a strange way



	all(types.>0) || error("Some input not understood")


	out1 = merge(hcat, to_myODict(), to_myODict.(x[types.==1], D)...) 


	function good_keys(used_keys, n, sol=[])
		
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


	return merge!(hcat, out1, 
								(OrderedDict(k=>to_myMatrix(v, D)) 
								 					for (k,v) in zip(unused_keys,x[types.==2]))...)


end





#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#




function to_myODict(x, arg::Union{<:Lattice,<:AbstractMatrix}
									 )::OrderedDict{Any, AbstractMatrix{Float64}}

	to_myODict(x, VecDim(arg))

end 


#to_myODict(x...)::OrderedDict{Any, AbstractMatrix{Float64}} = to_myODict(x)

#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

LattDims(latt::Lattice)::Vector{Int} = latt.LattDims

LattDims(latt::Lattice, d::Int)::Vector{Int} = LattDims(latt, [d])

function LattDims(latt::Lattice, S::Symbol)::Vector{Int} 

	S==:full || error("only works with S=:full")

	return 1:LattDim(latt.LattVec)

end 



function LattDims(latt::Lattice, 
									d::AbstractVector{<:Int})::AbstractVector{<:Int}

	@assert issubset(d,	LattDims(latt, :full)) "invalid operation"

	return d

end 


function LattDims(f::Function, latt::Lattice)::Vector{<:Int}

	d = copy(LattDims(latt))

	return applicable(f, d) ? LattDims(latt, f(d)) : d

end 

function LattDims!(f::Function, latt::Lattice)::AbstractVector{<:Int}

	d = LattDims(latt)

	!applicable(f, d) && return d 

	return LattDims!(latt, f(d))
	
end 


function LattDims!(latt::Lattice, d::Union{<:Int,AbstractVector{<:Int}}
									)::AbstractVector{<:Int}

	latt.LattDims = LattDims(latt, d)

end 




function LattVec(f::Function, latt::Lattice)::Matrix

	v = copy(LattVec(latt)) 
	
	return applicable(f, v) ? f(v) : v

end 

LattVec(::Nothing)::Nothing = nothing

LattVec(v::AbstractMatrix, x::Nothing=nothing)::AbstractMatrix = v

function LattVec(v::AbstractMatrix, D::AbstractVector{Int})::AbstractMatrix

	selectdim(v, 2, D)

end 

LattVec(latt::Lattice)::AbstractMatrix = LattVec(latt.LattVec, LattDims(latt))


function LattVec!(latt::Lattice, v::AbstractMatrix)::AbstractMatrix 

	setindex!(LattVec(latt), v, :, :)

end 

function LattVec!(f::Function, latt::Lattice)::AbstractMatrix

	v = LattVec(latt) 

	return applicable(f, v) ? LattVec!(latt, f(v)) : v 	

end 




#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#




LattDim(A::AbstractMatrix, x::Nothing=nothing)::Int = size(A,2)

LattDim(A::AbstractMatrix, D::AbstractVector{Int})::Int =LattDim(LattVec(A,D))

LattDim(latt::Lattice)::Int = LattDim(LattVec(latt))





VecDim(v::AbstractVecOrMat) = size(v,1)

function VecDim(latt::Lattice)::Int

	s = VecDim(LattVec(latt))

	s>0 && return s 

	for (prop, D) in latt, (k,v) in D 
	
		isnothing(v) && continue
			
		s = VecDim(v) 

		s>0 && return s 

	end 

	return 0

end 






#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function CombsOfVecs(arg::Union{Lattice, AbstractMatrix}, coeff
										 )::Matrix

	Utils.CombsOfVecs(LattVec(arg), 
										to_myMatrix(coeff, LattDim(arg));
										dim=2)
	
end


function CombsOfVecs(t::Tuple{<:Union{Lattice, AbstractMatrix}, Any};
										 kwargs...)::Matrix

	CombsOfVecs(t...; kwargs...)

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


function sublatt_labels(latt::Lattice, kind::Symbol)::Vector

	sublatt_labels(latt; kind=kind)

end 

	
function sublatt_labels(latt::Lattice; kind::Symbol=:Atoms)::Vector

	kind!=:Both && return sublatt_labels(latt[kind])

	return union(sublatt_labels(latt; kind=:Atoms),
							 sublatt_labels(latt; kind=:Vacancies))

end 

function sublatt_labels(D::OrderedDict)::Vector 

	collect(keys(D))

end 

function sublatt_labels(latts::Utils.List, args...; kwargs...)::Vector 

	Utils.isList(latts, Lattice) || error()

	return mapreduce(vcat, latts) do latt 

			sublatt_labels(latt, args...; kwargs...)

		end |> unique 

end 

#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function sublattices_contain(latt::Lattice, inp::Nothing=nothing; kwargs...)
	
	sublattices_contain(latt,""; kwargs...)

end 


function sublattices_contain(latt::Lattice, inp::Real; kwargs...)

	sublattices_contain(latt, string(Int(trunc(inp))); kwargs...)

end 


function sublattices_contain(latt::Lattice, inp::AbstractString; kwargs...)

	sublattices_contain(latt, [inp]; kwargs...)

end


function sublattices_contain(latt::Lattice, inp::Utils.List; kwargs...)
	
	# find all sublattices which contain items of inp 

	all_sl = sublatt_labels(latt; kwargs...)

	sl = filter(s->any(i->occursin(i,s), inp), all_sl) 

	return isempty(sl) ? all_sl : sl 

	#isempty(sl) && error("Wrong sublattice labels!")

end 


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


EmptyPos(D::Int64)::Matrix{Float64} = zeros(Float64, D, 0)

EmptyPos(n::Nothing)::Matrix{Float64} = EmptyPos(0)

function EmptyPos(arg::Union{<:Lattice, <:AbstractMatrix})::Matrix{Float64}
	
	EmptyPos(VecDim(arg))

end


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function PosAtoms(latt::Lattice;
									labels_contain=nothing,
									label=nothing,
#									f=(sl,P)->P,
									kind=:Atoms,
									kwargs...
									)

	SL = if !isnothing(label)
		
					in(label, sublatt_labels(latt, kind)) ? [label] : []
					
					else 
		
						sublattices_contain(latt, labels_contain)

				end 

	return mapreduce(s-> latt[kind][s], hcat, SL, init=EmptyPos(latt))

end






#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function UnitCells(arg::Union{<:Lattice,<:AbstractMatrix}, stopstart...)

	Utils.CombsOfVecs10(LattVec(arg), stopstart...; dim=2)

end 





#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function OuterOp(f::Symbol)::Function
	
	function outer_op(args::Vararg{Any,2}; kwargs...)::Array

		getproperty(Algebra, f)(map(args) do X 

			isa(X,Lattice) ? PosAtoms(X; kwargs...) : X

		end...; dim=2)

	end 
end


OuterDiff = OuterOp(:OuterDiff)
OuterDist = OuterOp(:OuterDist)
FlatOuterDist = OuterOp(:FlatOuterDist)
FlatOuterDiff = OuterOp(:FlatOuterDiff)
FlatOuterSum = OuterOp(:FlatOuterSum)


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#



function Distances(latt::Lattice, nr_neighbors::Int; 
									 nr_uc::Int=nr_neighbors, kwargs...)

	Distances(latt, :all; nr_uc=nr_uc, kwargs...)[1:nr_neighbors]

end 


function Distances(latt::Lattice, nr_neighbors::Symbol=:all; 
									 nr_uc::Int=2, kwargs...)

	haskey(kwargs,:nr_neighbors) && @warn "Obsolete kwarg"

	AtomsUC = PosAtoms(latt; kwargs...)
	
	Ds = FlatOuterDist(AtomsUC[:, 1:1].-AtomsUC, UnitCells(latt, nr_uc))

	Unique!(Ds; sorted=true)

	deleteat!(Ds, Ds.<TOLERANCE)

	nr_neighbors==:all && return Ds 

	error("Not understood: 'nr_neighbors'=$nr_neighbors")

#	return nr_neighbors=="all" ? Ds : Ds[1:min(nr_neighbors,end)]
  
end 


function BondDirs(latt::Lattice, neighbor_index::Int64=1; 
												nr_uc::Int=neighbor_index, kwargs...)

	AtomsUC = PosAtoms(latt; kwargs...)

	Diffs = FlatOuterDiff(AtomsUC[:, 1:1].-AtomsUC, UnitCells(latt, nr_uc))

	Dists,Inds = Unique(sum(abs2, Diffs,dims=1)[:], sorted=true, inds=:all) 

	choose = findall(Dists .> TOLERANCE)[neighbor_index]

	return Diffs[:, Inds[choose]]

end 



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#



function SquareInt_LattCoeff((latt,n)::Tuple{<:Lattice, <:Any})::Matrix{Int}

	SquareInt_LattCoeff(latt, n)

end
														

function SquareInt_LattCoeff(latt::Lattice, n)::Matrix{Int}

	n = to_myMatrix(n, latt)

	tn = trunc.(n)

	if !approx(tn, n)

		error("n is not an integer/an array of integers!")

	end

	all(size(n).==LattDim(latt)) && return tn 

	length(n)==LattDim(latt) && return LA.Diagonal(tn[:])

	error("Incorrect size of the matrix 'n' provided")


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




#===========================================================================#
#
# Obtain the ucs in big UC by checking which ucs are inside a polygon
#
#---------------------------------------------------------------------------#

function ucsUC_Polygon(v, mM, polygon; asfunction=false)

	dim=2 

	is_inside = Geometry.PointInPolygon_wn(polygon, 
																				 order_vertices=true, dim=dim)


	uc_candidates = Utils.vectors_of_integers(2, 	selectdim(mM, dim, 2).+1,
																								selectdim(mM, dim, 1).-1;
																								dim=dim)


	out(cand) = selectdim(cand, dim, 
												mapslices(is_inside, cand, dims=[2,1][dim])[:]
												) #|> collect


	!asfunction && return out(uc_candidates)



	shifts = sortslices(Utils.CombsOfVecs10(polygon, 1; dim=dim),
											dims=dim, by=LA.norm)

	return (
					
		size(shifts,dim),
					
		function shiftout(i)

			s = selectdim(shifts, dim, i:i)
	
			v = TOLERANCE * s / (LA.norm(s) + TOLERANCE/100)
	
			return out(v.+uc_candidates)

		end )



end 



#===========================================================================#
#
# Calculates the positons of the small unit cells inside the large Unit Cell
#		defined by the integer vectors in N
#
#---------------------------------------------------------------------------#

ucs_in_UC(N::Real; kwargs...) = ucs_in_UC(hcat(N); kwargs)

ucs_in_UC(N::Utils.List; kw...) = ucs_in_UC(LA.Diagonal(vcat(N...)); kw...) 

function ucs_in_UC(N::AbstractMatrix{T};  method2D_cells=:Polygon, kwargs...) where T

	if T<:Int 

		# Nothing to do, N is ready.

	elseif T<:Float64 || any(n->isa(n,Float64), N)

		tN = trunc.(N) 

		approx(tN, N) || error("N should contain integers")
		
		return ucs_in_UC(convert(AbstractMatrix{Int}, tN); kwargs...)

	else 

		error("Type '$T' not supported")

	end 


	correct_nr = Int(round(abs(LA.det(N))))
	# this should be the size of the new unit cell -- like | a1.(a2 x a3) |

	correct_nr>0 || error("Vectors must be linear independent")

	polygon = Geometry.rhomboidVertices_fromDVectors(N; dim=2)

	mM = cat(cat.(extrema(polygon, dims=2)..., dims=1)..., dims=2)


#  if len(N) == 1:
#
#    iter_ucs = [np.arange(*maxmin[::-1]).reshape(-1,1)]
#



	n, iter_ucs = if size(N,2)==2 && method2D_cells == :Polygon
									
									ucsUC_Polygon(N, mM, polygon; asfunction=true)

								else 

									error("Not implemented")

								end 

#    elif method2D_cells == "Spiral":
#      iter_ucs = ucsUC_Spiral(N,maxmin,correct_nr)
#
#    else:
#      raise ValueError("The second argument, 'method2D_cells', must be either 'Polygon' or 'Spiral'.")
#
#
#    
#  else:
#
#    iter_ucs = method2(N,maxmin)
#
#
	for i=1:n 

		ucs = iter_ucs(i)

		size(ucs, 2)==correct_nr && return ucs
		
		# if the correct number of uc-s in the large UC is found

	end 

	error("Could not construct big unit cell correctly.")

end

#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function parse_input_RsNs( latt::Union{Nothing,Lattice}=nothing;
													 Rs=nothing,
													 Ns=nothing,
													 A=nothing,
													 StopStart=nothing,
													 kwargs...
													)

	haskey(kwargs, :dim) && @warn "Obsolete kwarg dim"


	!isnothing(Rs) && return Rs 

	A = Utils.Assign_Value(A, LattVec(latt))

	isnothing(A) && error("Cannot proceed with unknown vectors")

	!isnothing(Ns) && return CombsOfVecs(A, to_myMatrix(Ns, A))

	!isnothing(StopStart) && return UnitCells(A, StopStart...)

	error("Not enough non-nothing kwargs")

end




#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

const labelCompStart = "_Part" 
const labelCompStop = "End" 

function labelToComponent(L::String)::String

	I,J = findall(labelCompStart,L), findall(labelCompStop,L)

	comp = Utils.mapif(!isnothing,  Base.product(I,J)) do (i,j) 

		a,b = i[end]+1, j[1]-1

		a<=b || return nothing 

		for other in (I,J)

			any(q->a<=q[1]<=b, other) && return nothing 

		end 

		return L[a:b] 

	end 



	if length(comp)>1 
		
		error("The name label '$L' contains too many '$labelCompStart...$labelCompStop' pairs")

	elseif length(comp)<1

		error("No component cue could be found in label '$L'")

	else 

		return comp[1]

	end 

end

function componentToLabel(L)::String

	string(labelCompStart,L,labelCompStop)

end 

function Labels_ManyUCs(P::Pair, uc_labels::AbstractVector)::Vector{String}

	Labels_ManyUCs(repeat([P.first], size(P.second,2)), uc_labels)

end

function Labels_ManyUCs(atom_labels::AbstractVector,
												uc_labels::AbstractVector;
#												kwargs...
												)::Vector{String}

	Algebra.OuterBinary(atom_labels, 
											componentToLabel.(uc_labels),
											*, flat=true)
end 



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#



function Atoms_ManyUCs(atoms::AbstractMatrix{<:Real},
											 ucs::AbstractMatrix{<:Real};
											 kwargs...)::Matrix{Float64}

	haskey(kwargs, :dim) && @warn "Obsolete kwarg dim"

	FlatOuterSum(atoms, ucs) 

end 



function Atoms_ManyUCs(atoms::AbstractMatrix{<:Real},
											 ns::AbstractMatrix{<:Real},
											 A::AbstractMatrix{<:Real};
											 kwargs...)::Matrix{Float64}

	haskey(kwargs, :dim) && @warn "Obsolete kwarg dim"
	Atoms_ManyUCs(atoms, parse_input_RsNs(;Ns=ns, A=A, kwargs...))

end 

function Atoms_ManyUCs(atoms::AbstractMatrix{<:Real},
											 ns::AbstractMatrix{<:Real},
											 latt::Lattice;
											 kwargs...)::Matrix{Float64}

	haskey(kwargs, :dim) && @warn "Obsolete kwarg dim"
	Atoms_ManyUCs(atoms, parse_input_RsNs(latt; Ns=ns, kwargs...))

end 





function Atoms_ManyUCs(latt::Lattice; kwargs...)::Matrix{Float64}

	haskey(kwargs, :dim) && @warn "Obsolete kwarg dim"
	Atoms_ManyUCs(PosAtoms(latt; kwargs...), parse_input_RsNs(latt; kwargs...))

end







#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function new_atoms_dict(latt::Lattice,
															 ns_ucs::AbstractMatrix,
															 Labels::Nothing=nothing)::Function

	function out(atoms::AbstractMatrix, k)::OrderedDict
			
		OrderedDict(k=>Atoms_ManyUCs(atoms, ns_ucs, latt)) 

	end 

end 



function new_atoms_dict(latt::Lattice,
												ns_ucs::AbstractMatrix,
												Labels::Function)::Function 

	labels_ucs = Labels.(eachcol(Int.(ns_ucs))) 

	return function out(atoms::AbstractMatrix, k)::OrderedDict
		
		labels = Labels_ManyUCs(k=>atoms, labels_ucs)

		return to_myODict((labels, Atoms_ManyUCs(atoms, ns_ucs, latt)), latt)

	end 

end 





#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#



function Superlattice!(latt::Lattice, n; Labels=nothing, kwargs...)::Lattice
	
	n = SquareInt_LattCoeff(latt, n)
	
	new_atoms = new_atoms_dict(latt, ucs_in_UC(n; kwargs...), Labels)

	for (kind,D) in latt, k in sublatt_labels(D) 

		merge!(D, new_atoms(pop!(D,k), k)) # may contain new labels != k
	
	end 

	LattVec!(latt, CombsOfVecs(latt, n))

	return latt 

end









function Superlattice(latt::Lattice, n; kw...)::Lattice

	Superlattice_(latt::Lattice, SquareInt_LattCoeff(latt, n); kw...)

end 
										


function Superlattice_(latt::Lattice, n; Labels=nothing, kw...)::Lattice
	
	Lattice(latt;
					act_on_vectors=(v::AbstractMatrix)->CombsOfVecs(v, n),
					act_on_atoms=new_atoms_dict(latt, ucs_in_UC(n; kw...), Labels),
				 )

end







#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#







function Superlattice(Components::Utils.List, Ns::Utils.List; Labels=nothing, kwargs...)::Lattice

	Ns = SquareInt_LattCoeff.(zip(Components,Ns))

	Supervectors = CombsOfVecs(Components[1], Ns[1])

	for S in CombsOfVecs.(zip(Components[2:end], Ns[2:end]))

		approx(Supervectors, S) && continue

		error("The 'Ns' provided are not compatible; one should recover the same supervectors for each pair (N,Component).")
	
	end



	Ns = ucs_in_UC.(Ns; kwargs...)

	Labels = 	if isnothing(Labels)
	
							repeat([""], length(Components))

						else 

							Utils.Assign_Value.(Labels, "")

						end 

	empty_labels = isempty.(Labels)


	return Lattice(Supervectors, map([:Atoms, :Vacancies]) do p

		Utils.flatmap(sublatt_labels(Components; kind=p)) do k 

			all_atoms = map(zip(Components,Ns)) do (latt,ns)

				Atoms_ManyUCs(latt; Ns=ns, label=k, kind=p)
			
			end 

			good_atoms = .!isempty.(all_atoms)

			together_labels = good_atoms .&   empty_labels
			separate_labels = good_atoms .& .!empty_labels
			

			return vcat(if !any(together_labels) 
									
										[] 
									
									else 

										 k=>hcat(all_atoms[together_labels]...)

									end,

									map(findall(separate_labels)) do i 

										string(Labels[i], k) => all_atoms[i]

									end)

		end  # Utils.flatmap ends 

	end...) # map over [:Atoms, :Vacancies] ends

end








#===========================================================================#
#
# 			some safety measureas when adding atoms 
#
#
#---------------------------------------------------------------------------#

function Initialize_PosAtoms(latt::Lattice, pos_atoms...)::AbstractMatrix

	Initialize_PosAtoms(VecDim(latt), pos_atoms...)

end 

function Initialize_PosAtoms(vect_dim::Int, pos_atoms=nothing)::AbstractMatrix
	
#    if type(pos_atoms)==type(None) or np.size(pos_atoms)==0:
#      print("\n *** No position given. The atom is set at the origin. ***\n")

	# if nothing is given, take the origin

	pos_atoms = to_myMatrix(if isnothing(pos_atoms) || isempty(pos_atoms) 

														zeros(Float64, vect_dim, 1)

													else 

														pos_atoms

													end)

	dim,nr = size(pos_atoms)

	dim==vect_dim && return pos_atoms
	
	dim>vect_dim && return vcat(pos_atoms, zeros(Float64, vect_dim-dim, nr))
	

	# in case additional coordinates are needed, use zeross

	error("Atoms cannot have more coordinates than the dimension of the lattice vectors!")

end 









#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function parse_input_sublattices(latt::Lattice, inp::Nothing, args...)

	good = !in(sublatt_labels(latt, :Both))

	for c in string.('A':'Z')

		good(c) && return parse_input_sublattices(latt, [c], args...)

	end 

	error()

end 

function parse_input_sublattices(latt::Lattice, inp::T, args...) where T<:Union{AbstractString, Char, Symbol, Real}

	parse_input_sublattices(latt, [inp], args...)

end 


function parse_input_sublattices(latt::Lattice, inp::Utils.List, nr_at, method)

	repeat(map(method, inp), outer=max(1,Int(ceil(nr_at/length(inp)))))[1:nr_at]
		# repeat the list if necessary

end 


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#



function AddAtoms!(latt::Lattice, positions=[], labels="A"; kind=:Atoms)

	positions = Initialize_PosAtoms(latt, positions)

	new_labels = parse_input_sublattices(latt, labels, size(positions,2), string)

	for (k,i) in zip(Unique(new_labels, inds=:all)...)

		isempty(i) && continue

		latt[kind][k] = hcat(get(latt[kind], k, EmptyPos(latt)), 
												 view(positions, :, i))

	end 

	return latt 

end 



function AddAtoms(latt::Lattice, positions=[], labels="A"; kind=:Atoms)

	positions = Initialize_PosAtoms(latt, positions)

	new_labels, inds = Unique(
												parse_input_sublattices(latt, labels, 
																								size(positions,2), string),
												inds=:all) 


	function act_on_atoms(K::Symbol, D::AbstractDict)
	
		kind!=K && return D
	
		return map(findall(!isempty, inds)) do not_empty

			k,R = new_labels[not_empty], view(positions, :, inds[not_empty])

			return Pair(k, haskey(D, k) ? hcat(D[k], R) : R)

		end 

	end 

	return Lattice(latt, act_on_atoms=act_on_atoms)

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

function parse_input_shift(latt::Lattice, R=nothing; 
													 n=nothing, r=nothing, kind=:Both)

	shift = if !isnothing(R) 
		
		to_myMatrix(latt, R)

					elseif !isnothing(r) 
					
		to_myMatrix(latt, r)

					else 

		parse_input_RsNs(latt; Ns=n, Rs=r)

					end 

	@assert size(shift,2)==1

	return function shift_atoms(K::Symbol, atoms::AbstractMatrix) 

		in(kind,[:Both,K]) ? atoms .+ shift : atoms 

	end 

end


function ShiftAtoms!(latt::Lattice, args...; kwargs...,)

	applyOnLattAtoms!(parse_input_shift(latt, args...; kwargs...), latt)
										
end


function ShiftAtoms(latt::Lattice, args...; kwargs...)

	Lattice(latt, act_on_atoms=parse_input_shift(latt, args...; kwargs...))

end 




#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

#function ReduceDim(latt::Lattice, dim=nothing; complement=false)
	
#	if complement && (isnothing(dim) || length(dim)>1)
#		
#		error("The complement can be computed only if one single dimension is given.")
#
#	end 
#
#	complement && return Lattice(latt, act_on_vectors=v->v[:, dim[1]:dim[1]])
#
#	isnothing(dim) && return Lattice(latt, act_on_vectors=EmptyPos)
#
#	keep_dims = filter(!in(vcat(dim...)),1:LattDim(latt))
#
#	return Lattice(latt, act_on_vectors=v->v[:, keep_dims])

function KeepDim!(latt::Lattice, args...)::Lattice 

	LattDims!(latt, args...)

	return latt 

end 


function KeepDim(latt::Lattice, args...)::Lattice 

	d = LattDims(latt, args...)

	return Lattice(latt; act_on_vectors=(::Vector{Int})->d)

end


function ReduceDim!(latt::Lattice, d=LattDims(latt))::Lattice 

	LattDims!(latt, setdiff(LattDims(latt), d))

	return latt

end 

function ReduceDim(latt::Lattice, d=LattDims(latt))::Lattice 

	LD = LattDims(latt, setdiff(LattDims(latt), d)) 

	return Lattice(latt; act_on_vectors=(::Vector{Int})->LD)

end 





#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


#function Distances(latt::Lattice; nr_uc=2, nr_neighbors=1, with_zero=true, kwargs...)
#
#	AtomsUC = PosAtoms(latt; kwargs...)
#	
#	Ds = FlatOuterDist(AtomsUC[:, 1:1].-AtomsUC, UnitCells(latt, nr_uc))
#
#	Utils.Unique!(Ds; tol=TOLERANCE, sorted=true)
#
#	!with_zero && deleteat!(Ds, Ds.<TOLERANCE)
#
#	return nr_neighbors=="all" ? Ds : Ds[1:min(nr_neighbors+1,end)]
#  
#end 

function get_Bonds(latt::Lattice; nr_uc=1, neighbor_index=1, kwargs...)::Vector{Tuple{Int,Int}}

	AtomsUC = PosAtoms(latt; kwargs...)

	size(AtomsUC,2) > 5000 && error("Not enough memory")

	UCs = UnitCells(latt, nr_uc)

	outside = map(!approx(zeros(VecDim(latt))), eachcol(UCs))

#	dists = zeros(size(AtomsUC,2), size(AtomsUC,2), size(UCs,2))
#
#	for (i,U) in enumerate(eachcol(UCs))
#
#		dists[:,:,i] = OuterDist(AtomsUC, AtomsUC .+ U)
#	
#	end 

	ud,ui = Unique(FlatOuterDist(AtomsUC, FlatOuterSum(AtomsUC,UCs)),
								 sorted=true, inds=:all)

#	@show approx(dists[:],FlatOuterDist(AtomsUC, FlatOuterSum(AtomsUC,UCs)))

	n = ui[neighbor_index + count(ud.<TOLERANCE)]

	ci = CartesianIndices((axes(AtomsUC,2),axes(AtomsUC,2),axes(UCs,2)))[n] 

	return Utils.mapif(!isempty, ci) do cart_idx 

		(i,j,u) = cart_idx.I 

		return outside[u] | isless(i,j) ? (i,j) : ()

	end |> unique |> sort 


end 





#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

#function SurfaceAtoms(latt::Lattice, bondlength=nothing; kwargs...)
#
#	LattDim(Latt)==0 || error()
#
#	SurfaceAtoms(PosAtoms(latt; kwargs...), bondlength)
#
#end 


#function SurfaceAtoms(atoms::AbstractMatrix{Float64})

#	size(atoms,2) > 5000 && error("Too many atoms")
#
#
#	Ds = OuterDist(atoms, atoms)
#
#	Utils.Uniqu1
#
#	UD = Utils.Unique(Ds; tol=TOLERANCE, sorted=true)[2:end]
#
#
#	bonds = Algebra.get_Bonds(atoms;
#														inds=true, pos=false, asmatrices=true)
#
#
#	bonds[
#
#

#convex hull


#end 





#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function Maximize_ContactSurface(starting_points, direction; vertices=nothing, ordered_vertices=nothing)
	

#    from shapely.geometry import LinearRing,MultiPoint,LineString
#    from shapely.ops import nearest_points
#  
#    if ordered_vertices is not None:
#        V = ordered_vertices
#    else:
#        V = Order_PolygonVertices(vertices)
#   
#    poly_sides,poly_verts = LinearRing(V),MultiPoint(V)
#  
#
#    shifts = np.zeros((2,2))
#
#    for start in starting_points:
#
#        line = LineString([start,start+direction])
#       
#        if poly_sides.intersection(line).is_empty:
#        
#            Ps = [p.coords for p in nearest_points(line,poly_verts)]
#
#            s = np.diff(Ps,axis=0).reshape(-1)
#            
#            for (i,shift) in enumerate(shifts):
#
#                if np.dot(shift,s) >= np.dot(shift,shift):
#
#                    shifts[i] = s
#
#    norms = np.linalg.norm(shifts,axis=1)
#
#    if any(norms<1e-10):
#        return shifts[np.argmax(norms),:]
#
#    return shifts[np.argmin(norms),:]
#
end 

function Align_toAtoms(latt::Lattice, atoms::AbstractMatrix{Float64}, shift_dir::Int=+1; bonds=nothing)

	# shift_dir = +/- i means positive/negative a_i
	# slowed down if all atoms provided instead of the surface ones

	if isnothing(bonds)
		
		c = Utils.DistributeBallsToBoxes.([-1,1], LattDim(latt))

		bonds = CombsOfVecs(latt, vcat(c...))

	end 

	atoms = Order_PolygonVertices(atoms)


#PointInPolygon_wn(SurfaceAtoms; order_vertices=false)




#    def is_inside(l):
#
#        for a in l.PosAtoms():
#    #            if Geometry.PointInPolygon_wn(a,SurfaceAtoms):
#
#                return True
#
#        return False
#
#
#    while is_inside(self):
#
#        self.Shift_Atoms(n=np.sign(-shift_dir),method="self")
#
#    while self.Dist_to_is(SurfaceAtoms,0).any():
#
#        self.Shift_Atoms(n=np.sign(-shift_dir),method="self")



	dmax = maximum(Algebra.OuterDist(atoms, atoms))*5


	dmax=10

	direction = sign(shift_dir)*LattVec(latt)[:,abs(shift_dir)]


	longdir = 2*dmax * LA.normalize(direction)


#	ShiftAtoms!(latt, r=-longdir/2)


	Maximize_ContactSurface(PosAtoms(latt),
													longdir,
													ordered_vertices=atoms)

#
#    self.Shift_Atoms(r=Geometry.Maximize_ContactSurface(
#                            self.PosAtoms(),
#                            longdir,
#                            ordered_vertices=SurfaceAtoms),
#                    method="self")


	return latt 

#
#    contact = Geometry.ClosestContact_LinesPolygon(
#		self.PosAtoms(),
#		longdir,
#                ordered_vertices=SurfaceAtoms)
#
#    rij = contact["stop"] - contact["start"]
#
#    out = self.Shift_Atoms(r=rij)
#
#	# Now the lead head (atom 'i') coincides with one surface atom 'j'
#	# We want to go back such that there i and j don't overlap anymore
#	# This will enable a bond between i and j
#	# A bond b is 'going back' if cos(b,rij) is as negative as possible
#
#    bonds = bonds[np.argsort(bonds/la.norm(bonds,axis=1,keepdims=True)@rij)]
#
#
#    def good(L):
#
#      for n in range(2):
#
#        if L.Dist_to_is(n*direction+SurfaceAtoms,0).any():
#
#          return False
#
#      return True
#
#
#
#    while not good(out):
#
#      for s in bonds:
# 
#        out2 = out.Shift_Atoms(r=s)
#
#        if good(out2):
#
#          return out2
#
#      out = out.Shift_Atoms(r=bonds[0])
#
#    return out
#


	nothing


end 



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

#function Wigner_Seitz(A::AbstractMatrix)
#
#	B = view(A, 1:minimum(size(A)),1:minimum(size(A)))
#
#	for n in 1:5
#
#		Vor = Voronoi(CombsOfVecs10(B, n))
#
#		regs = [Vor.vertices[r] for r in Vor.regions if !isempty(r)&&!in(-1,r)]
#
#		if !isempty(regs)
#  
##			reg = regs[argmin(map(LA.norm ∘ Algebra.Mean, regs))]
##VoronoiCells 
##			reg = regs[np.argmin(la.norm(np.mean(regs,axis=1),axis=1))]
#
#      return Geometry.Order_PolygonVertices(reg; dim=2)
#
#  raise Exception("The Wigner-Seitz cell couldn't be computed!")
#
#
#
#end 

#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function BrillouinZone(latt::Lattice)::Matrix

	LattDim(latt)==0 && return zeros(1,1)

  K = self.ReciprocalVectors() 

	LattDim(latt)==1 && return hcat(zero(K), K) 

	LattDim(latt)==2 && return Wigner_Seitz(K)

	error()

end 

#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function NearbyUCs(latt::Lattice, nr_uc::Int=1; #sublatt_coord=Dict(), 
									 remove_minus_m=true, kwargs...)::NTuple{3,Matrix}

	ms = Utils.vectors_of_integers(LattDim(latt), nr_uc, dim=2)

	inds = map(enumerate(eachcol(ms))) do (i,v)

		!remove_minus_m && return true 

		v==-v && return true 

		return !in(-v,eachcol(ms[:,1:i]))

	end 



	return (ms[:,inds],
					CombsOfVecs(latt, ms[:,inds]),
					PosAtoms(latt; kwargs...))






end 




#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function plot(latt::Lattice)

	Rs = PosAtoms(latt)

	Plots.scatter(Rs[1,:],Rs[2,:])

end 

function plot!(latt::Lattice;kwargs...)

	Rs = PosAtoms(latt)

	Plots.scatter!(Rs[1,:],Rs[2,:];kwargs...)

end 


function plot(latts...)

	plot(latts[1])

	for latt in latts[2:end-1]

		plot!(latt)

	end 

	atoms = hcat([PosAtoms(l) for l in latts]...)

	@show size(atoms)
#import 2libs 

	A = collect(eachcol(atoms))
@time h = ConcaveHull.concave_hull(A)
	h = hcat(h.vertices...)

	xlim,ylim = map(extrema.(eachrow(hcat(atoms,h)))) do lim 
		lim .+ (lim[2]-lim[1])*0.05*[-1,1]
								end

#	h = Order_PolygonVertices(h; dim=2)


B = collect(transpose(atoms))

@time CH = QHull.chull(B)


CH = transpose(B[CH.vertices,:])

#	Plots.plot!(h[1,:],h[2,:])
	Plots.plot!(CH[1,:],CH[2,:])




	plot!(latts[end],xlims=xlim,ylims=ylim,ratio=1)


#	extrema.(PosAtoms.(latts),dims=2)




#	return map(latts) do latt
#
#		Rs = PosAtoms(latt)
#
#		return Plots.scatter!(Rs[1,:],Rs[2,:])
#
#	end 
#


end 




#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function SquareLattice(name="A")

	Lattice(ArrayOps.UnitMatrix(2), name=>0)

end 




#############################################################################
end
