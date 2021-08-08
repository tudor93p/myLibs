module Lattices
#############################################################################

using OrderedCollections:OrderedDict

import ..LA

import ..Utils, ..Algebra, ..ArrayOps, ..Geometry

#import Plots, ConcaveHull, QHull



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


Base.copy(s::Symbol)::Symbol = s
Base.copy(s::AbstractString)::String = s




#===========================================================================#
#
# Functions that compute the 	-- number of vectors, 
# 														-- their indices 
# 														-- or get the vectors themselves
#
#  Lattice object not used
#
#---------------------------------------------------------------------------#

include("Lattices_VecMatDict.jl")

#===========================================================================#
#
# functions to_myMatrix And to_myODict 
#  Lattice object not used
#
#---------------------------------------------------------------------------#

# include("Lattices_MatDict.jl") # moved into the other file

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


		nr_vec, vec_len = NrVecs(lv, LD), VecLen(lv) 




		for k in [keys(sl);keys(va)], odict in [sl,va]

			if haskey(odict, k)

				s = VecLen(odict[k])

				if mode==:fractional

					s==nr_vec || error("Dimensionality mismatch. $s coefficient(s) provided for sublattice '$k', but $nr_vec are necessary")

					odict[k] = CombsOfVecs(Vecs(lv, LD), odict[k])
				
				elseif mode==:cartesian

					vec_len!=s && error("Dimensionality mismatch. The positions for sublattice '$k' are $s-dimensional, but the lattice vectors are $vec_len-dimensional")

				else 

					error("Mode $mode not defined")

				end 


			end 

		end 



		return new(lv, Utils.Assign_Value(LD, 1:nr_vec), sl, va)

	end 


end 


#===========================================================================#
#
# Simple square lattice 
#
#---------------------------------------------------------------------------#


function SquareLattice(name="A")::Lattice

	Lattice(ArrayOps.UnitMatrix(2), name=>0)

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
# iterate Lattice object
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
# Extract properties of Lattice object
#
#---------------------------------------------------------------------------#


function isintact(latt::Lattice)::Bool 

	NrVecs(latt.LattVec)==length(latt.LattDims)

end 

function dimsordered(latt::Lattice)::Bool 

	for (i,d) in enumerate(latt.LattDims) 

		i!=d && return false 

	end 

	return true  

end 



#===========================================================================#
#
# get Lattice dimension 
#
#---------------------------------------------------------------------------#

LattDim(latt::Lattice)::Int = length(latt.LattDims) 

function LattDim(latt::Lattice, mode::Symbol)::Int

	mode in (:free, :absolute, :relative) && return LattDim(latt)

	(mode==:stored || isintact(latt)) && return NrVecs(latt.LattVec)
	
	error("'mode=$mode' not supported")

end  



#===========================================================================#
#
# get Lattice vectors
#
#---------------------------------------------------------------------------#

LattVec(A::AbstractMatrix{<:Number})::Matrix{Number} = A 

LattVec(latt::Lattice)::Matrix{Float64} = Vecs(latt.LattVec, latt.LattDims) 

function LattVec(latt::Lattice, mode::Symbol)::Matrix{Float64}

	mode in (:free, :absolute, :relative) && return LattVec(latt) # DEFAULT 

	(mode==:stored || isintact(latt)) && return copy(latt.LattVec)

	error("'mode=$mode' not supported")

end 


#LattVec(::Nothing)::Nothing = nothing
#
#LattVec(v::AbstractMatrix, x::Nothing=nothing)::AbstractMatrix = v
#
#function LattVec(v::AbstractMatrix, D::AbstractVector{Int})::AbstractMatrix
#
#	Vecs(v, D) 
#
#end 


function LattVec(f::Function, latt::Lattice, args...)::Matrix{Float64}

	v = LattVec(latt, args...)
	
	return applicable(f, v) ? f(v) : v

end 




#===========================================================================#
#
# Change the lattice vectors and return the result 
#
#---------------------------------------------------------------------------#



function LattVec!(latt::Lattice, v::AbstractMatrix{<:Number}, 
									mode::Symbol=:free)::Matrix{Float64} 

	if mode==:stored || isintact(latt)  
		
		latt.LattVec = to_myMatrix(v)

	elseif mode in (:free, :absolute, :relative) 
		
		changed_inds = latt.LattDims 

		not_changed_inds = setdiff(IndsVecs(latt.LattVec), changed_inds)

		if maximum(changed_inds) < minimum(not_changed_inds)  

			latt.LattVec = Cat(v, LattVec(latt, not_changed_inds, :stored))

		elseif maximum(not_changed_inds) < minimum(changed_inds)  

			latt.LattVec = Cat(LattVec(latt, not_changed_inds, :stored), v)

		elseif NrVecs(v)==length(changed_inds)

			SetVecs!(latt.LattVec, changed_inds, v) 

		else 

			error("Cannot change only vectors $changed_inds")

		end 

	else 

		error("'mode=$mode' not supported")

	end 

	return copy(latt.LattVec)

end  



function LattVec!(f::Function, latt::Lattice, args...)::Matrix{Float64}

	v = LattVec(latt, args...)

	return applicable(f, v) ? LattVec!(latt, f(v), args...) : v 	

end 



#===========================================================================#
#
# get Lattice dimensions
#
#---------------------------------------------------------------------------#

LattDims(latt::Lattice)::Vector{Int} = 1:LattDim(latt) 


function LattDims(latt::Lattice, mode::Symbol)::Vector{Int} 

	mode==:relative && return LattDims(latt) # DEFAULT

	(mode==:stored || isintact(latt)) && return IndsVecs(latt.LattVec)

	(mode==:absolute || dimsordered(latt)) && return copy(latt.LattDims)

	error("'mode=$mode' not understood")

end 



# check if indices make sense and return them back 
 
function LattDims(latt::Lattice, d::Union{Int,AbstractVector{Int}}, 
									args...; complement::Bool=false)::Vector{Int} 

	D = LattDims(latt, args...)

	d_ = IndsVecs(D, d)

	return complement ? setdiff(D, d_) : d_ 

end 





function LattDims(f::Function, latt::Lattice, args...; 
									complement::Bool=false)::Vector{Int}

	D = LattDims(latt, args...)

	d_ = applicable(f, D) ? IndsVecs(D, f(D))  : D 

	return complement ? setdiff(D, d_) : d_ 

end 





#===========================================================================#
#
# Helper function for dim change
#
#---------------------------------------------------------------------------#

function is_mode_rel(args...)::Bool

	i_mode = findfirst(a->isa(a,Symbol), args)

	i_latt = findfirst(a->isa(a,Lattice), args)


	isnothing(i_mode) && return true   #not found => default==relative

	mode = args[i_mode]

	mode==:relative && return true


	mode in (:stored, :absolute) && return false 


	for f in (isintact, dimsordered)
	
		!isnothing(i_latt) && f(args[i_latt]) && return false 
	
	end 

	error("'mode=$mode' not understood")

end 


function new_LattDims(ld::Vector{Int}, latt::Lattice, args...)::Vector{Int}
	
	is_mode_rel(latt, args...) ? latt.LattDims[ld] : ld 

end 

#===========================================================================#
#
# Change latt.LattDims to d and return d
#
#---------------------------------------------------------------------------#

function LattDims!(latt::Lattice, args...; kwargs...)::Vector{Int} 

	ld = LattDims(latt, args...; kwargs...)

	latt.LattDims = new_LattDims(ld, latt, args...)

end 

function LattDims!(f::Function, latt::Lattice, args...; kwargs...
									 )::Vector{Int}

	hasmethod(f, (AbstractVector{Int},)) || return latt.LattDims  

	D = LattDims(latt, args...; kwargs...)
		
	latt.LattDims = new_LattDims(IndsVecs(D, f(D)), latt, args...)

end 



#===========================================================================#
#
# Extend the functions from "Lattices_Vecs.jl" and "Lattices_MatDict.jl"
#
#---------------------------------------------------------------------------#


#NrVecs 
#IndsVecs 
#Vecs 
#VecLen 
#VecAx 
#Shape 
#ZeroVecs 
#EmptyVec 
#


function VecLen(latt::Lattice)::Int

	s = VecLen(LattVec(latt, :stored))

	s>0 && return s 

	for (prop, D) in latt, (k,v) in D 
	
		isnothing(v) && continue
			
		s = VecLen(v) 

		s>0 && return s 

	end 

	return 0

end 


function EmptyVec(latt::Lattice)::Matrix{Float64}
	
	EmptyVec(VecLen(latt))

end




function ZeroVecs(latt::Lattice, args...)::Matrix{Float64}

	ZeroVecs(VecLen(latt), args...)

end 

function ZeroVecs(vec_len::Int, latt::Lattice, args...)::Matrix{Float64}

	ZeroVecs(vec_len, LattDim(latt, args...))

end 




function to_myMatrix(x, latt::Lattice)::Matrix{Float64}

	to_myMatrix(x, VecLen(latt))
	
end 


function to_myODict(x, latt::Lattice)::OrderedDict{Any, Matrix{Float64}}

	to_myODict(x, VecLen(latt))

end 


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#




function applyOnLattVec!(f::Function, latt::Lattice, args...)::Lattice 

	LattDims!(f, latt, args...)  # if applicable -- checked

	LattVec!(f, latt, args...)		# if applicable -- checked 

	# in LattDims,  f is applied on the dims of the actual latt.LattVec,
	# but in LattVec, f is applied only on the free vectors among latt.LattVec


	if applicable(f, latt) 

		v, d = f(latt)

		LattDims!(latt, d, args...)

		LattVec!(latt, v, args...)

	end 

	return latt 

end 



function applyOnLattVec(f::Function, latt::Lattice, args...
												)::Tuple{Matrix{Float64}, Vector{Int}}

	(LattVec(f, latt, args...), LattDims(f, latt, args...))

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

		applyOnLattAtoms(f, latt, kind, D)

	end
	
end 



function applyOnLattAtoms(f::Function, latt::Lattice, 
													kind::Symbol, D::AbstractDict=latt[kind])

	for args in ((D,), (kind,), (kind, D),) 

		applicable(f, args...) && return f(map(copy, args)...)

	end 

	return map(sublatt_labels(D)) do k 

		applyOnLattAtoms(f, latt, kind, k, D[k]) 

	end 

end 


function applyOnLattAtoms(f::Function, red_op::Function, latt::Lattice, 
													kind::Symbol, D::AbstractDict=latt[kind])

	mapreduce(red_op, sublatt_labels(D)) do k 

		applyOnLattAtoms(f, latt, kind, k, D[k]) 

	end 

end 


function applyOnLattAtoms(f::Function, latt::Lattice, 
													kind::Symbol, k::Any,
													atoms::AbstractMatrix=latt[kind][k])

	for args in ((atoms,), (kind, atoms), 
							 (atoms, k), (kind, atoms, k),
							 (kind, k),
							 )

		applicable(f, args...) && return f(map(copy, args)...)

	end 

	error("No method of 'f' could be used")

end 


#function applyOnLattAtoms(f::Function, #latt::Lattice, 
#													kind::Symbol, atoms::AbstractMatrix)
#	
#	for args in ((atoms,), (kind, atoms))
#
#		applicable(f, args...) && return f(map(copy, args)...)
#
#	end 
#
#	error("No method of 'f' could be used")
#
#end 
#





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



function Lattice(latt::Lattice, args...; 
								 act_on_vectors::Function=()->nothing,#identity,
								 act_on_atoms::Function=identity,
																)

	vectors,dims = applyOnLattVec(act_on_vectors, latt, args...)

	return Lattice(vectors, applyOnLattAtoms(act_on_atoms, latt)..., dims)

end 
	
		

#===========================================================================#
#
#
#
#---------------------------------------------------------------------------# 





#===========================================================================#
#
# Combination of Vectors in an array or the latt vect of a Lattice 
#
#---------------------------------------------------------------------------#


function CombsOfVecs(latt::Lattice, coeff, args...)::Matrix{Float64}

	CombsOfVecs(LattVec(latt, args...), coeff)

end 


function CombsOfVecs(t::Tuple; kwargs...)::Matrix{Float64}

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

	filter(sublatt_labels(latt; kwargs...)) do s
	
		any(inp) do i 

			occursin(string(i), string(labelToComponent(s))) 

		end 

	end 

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

function parse_input_labels(latt::Lattice, kind::Symbol, labels_contain, label::Nothing=nothing)::Vector

	sublattices_contain(latt, labels_contain, kind=kind)

end 


function parse_input_labels(latt::Lattice, kind::Symbol, labels_contain, label)::Vector

	in(label, sublatt_labels(latt, kind)) ? [label] : []

end 


PosAtoms(X::AbstractArray; kwargs...)::Array = X

function PosAtoms(latt::Lattice;
									labels_contain=nothing,
									label=nothing,
									kind=:Atoms,
									kwargs...)::Matrix

	SL = parse_input_labels(latt, kind, labels_contain, label) 

	E = EmptyVec(latt) 

	isempty(SL) && return E 

	return applyOnLattAtoms(hcat, latt, kind) do atoms::AbstractMatrix,k::Any

		in(k, SL) ? atoms : E 

	end 

end 




function labelToComponent(latt::Lattice;
													labels_contain=nothing,
													label=nothing,
													kind=:Atoms,
													kwargs...)::Vector{String}
	
	SL = parse_input_labels(latt, kind, labels_contain, label) 

	E = String[]

	isempty(SL) && return E

	applyOnLattAtoms(vcat, latt, kind) do atoms::AbstractMatrix,k::String

		in(k,SL) ? labelToComponent(k, atoms) : E

	end


end 





#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function UnitCells(arg::Union{<:Lattice,<:AbstractMatrix}, stopstart...)

	CombsOfVecs10(LattVec(arg), stopstart...)

end 





#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function OuterOp(f::Symbol)::Function

	F = getproperty(Algebra, f) 

	return function outerop(args::Vararg{<:AbstractArray,2}; kwargs...)::Array

		F(PosAtoms.(args; kwargs...)...; dim=VECTOR_STORE_DIM)

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



function SquareInt_LattCoeff(t::Tuple)::Matrix{Int}

	SquareInt_LattCoeff(t...)

end
														

function SquareInt_LattCoeff(latt::Lattice, n, args...)::Matrix{Int}

	n = to_int(to_myMatrix(n, latt))

	ld = LattDim(latt, args...)

	length(n)==ld && return LA.Diagonal(n[:]) 

	LA.checksquare(n)==ld && return n

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

	is_inside = PointInPolygon(polygon, order_vertices=true)

	uc_candidates = vectors_of_integers(2, Vecs(mM, 2)[:].+1, Vecs(mM, 1)[:].-1)


	function out(cand::AbstractMatrix{T})::Matrix{T} where T<:Number
		
		Vecs(cand, map(is_inside, eachvec(cand)))

	end 


	!asfunction && return out(uc_candidates)

	shifts = sortvecs(CombsOfVecs10(polygon, 1), by=LA.norm)


	return (
					
		NrVecs(shifts),
					
		function shiftout(i::Int)

			s = Vecs(shifts, i)
	
			v = TOLERANCE * s / (LA.norm(s) + TOLERANCE/100)
	
			return out( v .+ uc_candidates )

		end )



end 


function to_int(x::AbstractArray{Int,N})::Array{Int,N} where N 

	x

end 


function to_int(x::AbstractArray{<:Number,N})::Array{Int,N} where N 

	tx = trunc.(x)

	@assert approx(tx, x) "x should contain integers"

	return tx

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

		N[:,:] .= to_int(N)

	else 

		error("Type '$T' not supported")

	end 


	correct_nr = Int(round(abs(LA.det(N))))
	# this should be the size of the new unit cell -- like | a1.(a2 x a3) |

	correct_nr>0 || error("Vectors must be linear independent")

	polygon = Geometry.rhomboidVertices_fromDVectors(N; dim=VECTOR_STORE_DIM)

	mM = Cat(cat.(extrema(polygon, dims=VECTOR_STORE_DIM)..., dims=VECTOR_AXIS)...)


#  if len(N) == 1:
#
#    iter_ucs = [np.arange(*maxmin[::-1]).reshape(-1,1)]
#



			n, iter_ucs = if NrVecs(N)==2 && method2D_cells == :Polygon
									
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
#
#    iter_ucs = method2(N,maxmin)
#
#
	for i=1:n 

		ucs = iter_ucs(i)

		NrVecs(ucs)==correct_nr && return ucs
		
		# if the correct number of uc-s in the large UC is found

	end 

	error("Could not construct big unit cell correctly.")

end

#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function parse_input_RsNs( latt::Union{Nothing,Lattice}=nothing, args...;
													 Rs=nothing, Ns=nothing, A=nothing, 
													 StopStart=nothing, kwargs...
													)

	haskey(kwargs, :dim) && @warn "Obsolete kwarg dim"


	!isnothing(Rs) && return Rs 

	A = Utils.Assign_Value(A, LattVec(latt, args...))

	isnothing(A) && error("Cannot proceed with unknown vectors")

	!isnothing(Ns) && return CombsOfVecs(A, to_myMatrix(Ns, NrVecs(A)))

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

		return L
		error("No component cue could be found in label '$L'")

	else 

		return comp[1]

	end 

end


function labelToComponent(L::T)::T where T
	
	L

end 


function labelToComponent(L, atoms::AbstractMatrix)::Vector

	repeat([labelToComponent(L)], size(atoms,2))

end 



function componentToLabel(L)::String

	string(labelCompStart, L, labelCompStop)

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

	return FlatOuterSum(atoms, ucs) 

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
											 latt::Lattice, args...;
											 kwargs...)::Matrix{Float64}

	haskey(kwargs, :dim) && @warn "Obsolete kwarg dim"
	
	Atoms_ManyUCs(atoms, parse_input_RsNs(latt, args...; Ns=ns, kwargs...))

end 





function Atoms_ManyUCs(latt::Lattice, args...; kwargs...)::Matrix{Float64}

	haskey(kwargs, :dim) && @warn "Obsolete kwarg dim" 

	Atoms_ManyUCs(PosAtoms(latt; kwargs...), 
								parse_input_RsNs(latt, args...; kwargs...))

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

	labels_ucs = Labels.(eachvec(Int.(ns_ucs))) 

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



function Superlattice!(latt::Lattice, n, args...;
											 Labels=nothing, kwargs...)::Lattice

	n = SquareInt_LattCoeff(latt, n, args...)
	
	new_atoms = new_atoms_dict(latt, ucs_in_UC(n; kwargs...), Labels)

	for (kind,D) in latt, k in sublatt_labels(D) 

		merge!(D, new_atoms(pop!(D,k), k)) # may contain new labels != k
	
	end 

	LattVec!(latt, CombsOfVecs(latt, n), args...)

	return latt 

end









function Superlattice(latt::Lattice, n, args...; kw...)::Lattice

	Superlattice_(latt, SquareInt_LattCoeff(latt, n, args...), args...; kw...)

end 
										


function Superlattice_(latt::Lattice, n::Matrix{Int}, args...;
											 Labels=nothing, kw...)::Lattice

	cov = CombsOfVecs(latt, n, args...)

	return Lattice(latt, args...;
					act_on_vectors=(v::AbstractMatrix)->cov,
					act_on_atoms=new_atoms_dict(latt, ucs_in_UC(n; kw...), Labels),
				 )

end







#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#







function Superlattice(Ls::Utils.List, Ns::Utils.List, args...; 
											Labels=nothing, kwargs...)::Lattice 

	Labels_ = Dict{Symbol,Vector{String}}()

	if isnothing(Labels) 

		Labels_[:Labels] = repeat([""], length(Ls))

	elseif Utils.isList(Labels) && length(Labels)==length(Ls)
		
		Labels_[:Labels] = [isnothing(L) ? "" : string(L) for L in Labels]

	else 

		error("Labels not understood")

	end 

	Ns_ = (SquareInt_LattCoeff(L, N, args...) for (L,N) in zip(Ls, Ns))

	return Superlattice_(collect(Ls), collect(Ns_), 
											 Labels_[:Labels], args...; kwargs...)

end 




function Superlattice_(Components::AbstractVector{Lattice},
											 Ns::AbstractVector{<:AbstractMatrix{Int}},
											 Labels::AbstractVector{<:AbstractString},
											 args...; kwargs...)::Lattice


	Supervectors = CombsOfVecs(Components[1], Ns[1], args...)

	for (c,n) in zip(Components[2:end], Ns[2:end])

		approx(Supervectors, CombsOfVecs(c, n, args...)) && continue

		error("The 'Ns' provided are not compatible; one should recover the same supervectors for each pair (N,Component).")
	
	end



	Ns = ucs_in_UC.(Ns; kwargs...)


	empty_labels = isempty.(Labels)



	return Lattice(Supervectors, map([:Atoms, :Vacancies]) do p

		Utils.flatmap(sublatt_labels(Components; kind=p)) do k 

			all_atoms = map(zip(Components,Ns)) do (latt,ns)

				Atoms_ManyUCs(latt, args...; Ns=ns, label=k, kind=p)
			
			end 

			good_atoms = .!isempty.(all_atoms)

			together_labels = good_atoms .&   empty_labels
			separate_labels = good_atoms .& .!empty_labels
			
			
			pairs_ = [string(Labels[i], k)=>all_atoms[i] 
																for i in findall(separate_labels)]


			!any(together_labels) && return pairs_ 

			return [k=>Cat(all_atoms[together_labels]...); pairs_]


		end  # Utils.flatmap ends 

	end...) # map over [:Atoms, :Vacancies] ends

end








#===========================================================================#
#
# 			some safety measureas when adding atoms 
#
#---------------------------------------------------------------------------#

function Initialize_PosAtoms(latt::Lattice, pos_atoms...)::Matrix{Float64}

	Initialize_PosAtoms(VecLen(latt), pos_atoms...)

end 

function Initialize_PosAtoms(vect_dim::Int, x::Nothing=nothing
														)::Matrix{Float64}
	
#    if type(pos_atoms)==type(None) or np.size(pos_atoms)==0:
#      print("\n *** No position given. The atom is set at the origin. ***\n")

	# if nothing is given, take the origin


	Initialize_PosAtoms(vect_dim, ZeroVecs(vect_dim, 1))

end 

function Initialize_PosAtoms(vect_dim::Int,
														 pos_atoms::AbstractArray{<:Number}
														 )::Matrix{Float64}

	isempty(pos_atoms) && return Initialize_PosAtoms(vect_dim)

	pos_atoms = to_myMatrix(pos_atoms, vect_dim)

	dim = VecLen(pos_atoms)

	dim==vect_dim && return pos_atoms
	
	dim>vect_dim && return cat(pos_atoms, 
														 ZeroVecs(vect_dim-dim, NrVecs(pos_atoms)), 
														 dims=VECTOR_AXIS)
	
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

	for (k,i) in EnumUnique(new_labels)
		
		isempty(i) && continue

		latt[kind][k] = hcat(get(latt[kind], k, EmptyVec(latt)), 
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

	new_dict = map(findall(!isempty, inds)) do not_empty

			new_labels[not_empty]=>view(positions, :, inds[not_empty])

	end |> OrderedDict


	function act_on_atoms(K::Symbol, D::AbstractDict)
	
		kind==K ? merge(hcat, D, new_dict) : D
		
	end 

	return Lattice(latt, act_on_atoms=act_on_atoms)

end


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function RemoveAtoms!(latt::Lattice,
											inds::AbstractVector{Int},
											sublatt)::Lattice

	isempty(inds) && return latt
	
	latt[:Atoms][sublatt], latt[:Vacancies][sublatt] = [latt[:Atoms][sublatt][:, i] for i in [setdiff(axes(latt[:Atoms][sublatt],2), inds), inds]]

	return latt 

end 



function RemoveAtoms!(latt::Lattice, 
											inds::AbstractVector{Bool},
											sublatt)::Lattice

	any(inds) || return latt 

	latt[:Atoms][sublatt],latt[:Vacancies][sublatt] = [latt[:Atoms][sublatt][:,i] for i in [.!inds, inds]] # unstable , bounds error 

	return latt 

end 

function RemoveAtoms(latt::Lattice,
										 inds::AbstractVector{Bool},
										 sublatt)::Lattice

	function act_on_atoms(K::Symbol, D::AbstractDict)
	
		K==:Atoms || return D

		D[sublatt] = D[sublatt][:,.!inds]

		return D  

	end 

	return Lattice(latt, act_on_atoms=act_on_atoms)

end 




#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function parse_input_shift(latt::Lattice, R=nothing; 
													 n=nothing, r=nothing, kind=:Both)

	shift = if !isnothing(R) 
		
		to_myMatrix(R, latt)

					elseif !isnothing(r) 
					
		to_myMatrix(r, latt)

					else 

		parse_input_RsNs(latt; Ns=n, Rs=r)

					end 

	@assert NrVecs(shift)==1

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
#	isnothing(dim) && return Lattice(latt, act_on_vectors=EmptyVec)
#
#	keep_dims = filter(!in(vcat(dim...)),1:LattDim(latt))
#
#	return Lattice(latt, act_on_vectors=v->v[:, keep_dims])



function KeepDim!(latt::Lattice, args...)::Lattice 

	LattDims!(latt, args...)

	return latt 

end 

# 

function KeepDim(latt::Lattice, args...)::Lattice 

	act(::AbstractVector{Int}) = LattDims(latt, args...) 

	return Lattice(latt, :stored; act_on_vectors=act) 

end



#===========================================================================#
#
# Reduce the lattice to zero-dim 
#
#---------------------------------------------------------------------------#

function ReduceDim!(latt::Lattice)::Lattice 

	LattDims!(latt, Int[])

	return latt

end 


function ReduceDim(latt::Lattice)::Lattice 
	
	act(::AbstractVector{Int})::Vector{Int} = Int[] 

	return Lattice(latt, :stored; act_on_vectors=act)

end  



#===========================================================================#
#
# Custom reduce dim 
#
#---------------------------------------------------------------------------#


function ReduceDim!(latt::Lattice, d::Union{Int,AbstractVector{Int}}, 
										args...)::Lattice 

	ld = LattDims(latt, d, args...; complement=true)

	LattDims!(latt, ld, args...)

	return latt

end  




function ReduceDim(latt::Lattice, d::Union{Int,AbstractVector{Int}}, 
										args...)::Lattice 

	ld = LattDims(latt, d, args...; complement=true)

	act(::AbstractVector{Int})::Vector{Int} = new_LattDims(ld, latt, args...)
	
	return Lattice(latt, :stored; act_on_vectors=act) 

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

function get_Bonds(latt::Lattice; 
									 nr_uc::Int=1, 
									 neighbor_index::Int=1, kwargs...)::Vector{Tuple{Int,Int}}

	AtomsUC = PosAtoms(latt; kwargs...)

	size(AtomsUC,2) > 5000 && error("Not enough memory")

	UCs = UnitCells(latt, nr_uc)

	outside = map(!approx(ZeroVecs(latt, 1)), eachvec(UCs))


	ud,ui = Unique(FlatOuterDist(AtomsUC, FlatOuterSum(AtomsUC,UCs)),
								 sorted=true, inds=:all)

#	@show approx(dists[:],FlatOuterDist(AtomsUC, FlatOuterSum(AtomsUC,UCs)))

	n = ui[neighbor_index + count(ud.<TOLERANCE)]

	ci = CartesianIndices((IndsVecs(AtomsUC), 
												 IndsVecs(AtomsUC), 
												IndsVecs(UCs)))[n] 

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

function SurfaceAtoms(latt::Lattice, args...; kwargs...)::Matrix{Float64}

	LattDim(latt)==0 || error()

	return SurfaceAtoms(PosAtoms(latt; kwargs...), args...)

end 


function SurfaceAtoms(atoms::AbstractMatrix{Float64}, 
											bondlength::Real)::Matrix{Float64}

	NrVecs(atoms) > 5000 && error("Too many atoms")

#	Ds = OuterDist(atoms, atoms)

#	d_nn = Utils.Unique(Ds[:]; tol=TOLERANCE, sorted=true)[2]

	nr_bonds = NrBonds(atoms, bondlength)
	
	return Vecs(atoms, nr_bonds .< maximum(nr_bonds))
	
end 

function NrBonds(latt::Lattice, args...; kwargs...)::Vector{Int}

	NrBonds(PosAtoms(latt; kwargs...), args...)

end 

function NrBonds(atoms::AbstractMatrix{Float64}, bondlength::Real)::Vector{Int} 
	out = zeros(Int, NrVecs(atoms))
	
	NrVecs(atoms) > 1 || return out 

	allpairs = Algebra.get_Bonds_asMatrix(atoms, bondlength; 
																				dim=VECTOR_STORE_DIM, inds=true) 

	isempty(allpairs) && return out


	for (i,bonds) in Utils.EnumUnique(allpairs[:])
		
		out[i] = length(bonds)

	end 

	return out 

end 


#function NrBonds(atoms1::AbstractMatrix{Float64}, 
#								 atoms2::AbstractMatrix{Float64},
#								 bondlength::Real)::Tuple{Vector{Int},Vector{Int}}
#
#	EnumUnique(allpairs[:,1])
#
#	out1[i] = ...
#	allpairs[:,2]
#
#end  



function RemoveSingleBonds!(Latt::Lattice, d_nn::Real)::Lattice


#sublatt_labels(Latt)::Vector

#	x = applyOnLattAtoms(hcat, Latt, :Atoms) do atoms::AbstractMatrix,k::Any
#
#		@show k 
#		axes(atoms,2)
#
##		in(k, SL) ? atoms : E 
#
#	end 
#
#	@show x 

	rmv_inds = NrBonds(Latt, d_nn) .<= 1

	any(rmv_inds) || return Latt

	RemoveAtoms!(Latt, rmv_inds, "A") #unstable 

	return RemoveSingleBonds!(Latt, d_nn)

end 


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function do_overlap(latt::Lattice, args...;kwargs...)::Bool

	do_overlap(PosAtoms(latt), args...; kwargs...)

end 


function do_overlap(atoms::AbstractMatrix{<:Real}, 
										surface::AbstractMatrix{<:Real}; 
										kwargs...)::Bool

	for P in eachvec(atoms) 

		if PointInPolygon(P, surface; prepare_vertices=false)

			return true

		end 

	end 

	return any(approx.(OuterDist(surface, atoms), 0))

end 



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function Align_toAtoms!(latt::Lattice, latt2::Lattice, args...;
												kwargs...)::Lattice

	Align_toAtoms!(latt, PosAtoms(latt2), Distances(latt2)[1], args...;
								 kwargs...)

end 


function Align_toAtoms!(latt::Lattice, 
												atoms::AbstractMatrix{Float64},
												d_nn::Float64, args...;
												prepare_vertices=false, order_vertices=false,
												kwargs...)::Lattice

	Align_toAtoms!(latt::Lattice, SurfaceAtoms(atoms, d_nn), args...; 
								 prepare_vertices=true, order_vertices=true, kwargs...) 
								 
end 



function Align_toAtoms!(latt::Lattice, 
												surface::AbstractMatrix{Float64}, 
												shift_dir::Int=-1; 
												kwargs...)::Lattice 
	
	Align_toAtoms!(latt, 
								 surface, 
								 sign(shift_dir)*LattVec(latt)[:,abs(shift_dir)];
								 kwargs...)

end 


function Align_toAtoms!(latt::Lattice, 
												surface::AbstractMatrix{Float64}, 
												shift_dir::AbstractVector{Float64};
												prepare_vertices=true, order_vertices=true, 
												kwargs...)::Lattice

	# shift_dir = +/- i means positive/negative a_i
	# slowed down if all atoms provided instead of the surface ones

#	if isnothing(bonds)
#		
#		c = Utils.DistributeBallsToBoxes.([-1,1], LattDim(latt))
#
#		bonds = CombsOfVecs(latt, hcat(vcat(c...)...))
#
#	end 
	
	prepare_vertices && return Align_toAtoms!(latt, Geometry.prepare_polygon_vertices(surface; dim=VECTOR_STORE_DIM, order_vertices=order_vertices), shift_dir; prepare_vertices=false, kwargs...)




	while true 

		r = Algebra.LargestLengthscale(surface; dim=VECTOR_STORE_DIM)

		ShiftAtoms!(latt; r=-shift_dir*r*1.3)
		
		do_overlap(latt, surface) || break 

	end 


	ShiftAtoms!(latt, r=Geometry.Maximize_ContactSurface(PosAtoms(latt),
																											 shift_dir,
																											 surface, 
																											 dim=VECTOR_STORE_DIM,
																											 prepare_vertices=false)
							)

	#bridge may be needed !! bonds must be provided  



	return latt 

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


function ReciprocalVectors(latt::Lattice)::Matrix{Float64}

	LV = LattVec(latt, :stored)

	LA.checksquare(LV)

	return Vecs(2pi*inv(transpose(LV)), LattDims(latt, :absolute))

end 

	

#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function BrillouinZone(latt::Lattice)::Matrix{Float64}

	LattDim(latt)==0 && return zeros(0,1)

	K = ReciprocalVectors(latt)

	LattDim(latt)==1 && return VecsToMat(zero(K), K)
																									
	error()

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

	ms = vectors_of_integers(LattDim(latt), nr_uc)

	ms1 = Vecs(ms, map(enumerate(eachvec(ms))) do (i,v)

		remove_minus_m || return true 

		v==-v && return true 

		return !in(-v, eachvec(Vecs(ms, 1:i)) )

	end)


	return (ms1, CombsOfVecs(latt, ms1), PosAtoms(latt; kwargs...))

end 




#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


#function plot(latt::Lattice)
#
#	Rs = PosAtoms(latt)
#
#	Plots.scatter(Rs[1,:],Rs[2,:])
#
#end 
#
#function plot!(latt::Lattice;kwargs...)
#
#	Rs = PosAtoms(latt)
#
#	Plots.scatter!(Rs[1,:],Rs[2,:];kwargs...)
#
#end 
#
#
#function plot(latts...)
#
#	plot(latts[1])
#
#	for latt in latts[2:end-1]
#
#		plot!(latt)
#
#	end 
#
#	atoms = hcat([PosAtoms(l) for l in latts]...)
#
#	@show size(atoms)
##import 2libs 
#
#	A = collect(eachvec(atoms))
#@time h = ConcaveHull.concave_hull(A)
#	h = hcat(h.vertices...)
#
#	xlim,ylim = map(extrema.(eachrow(hcat(atoms,h)))) do lim 
#		lim .+ (lim[2]-lim[1])*0.05*[-1,1]
#								end
#
##	h = Order_PolygonVertices(h; dim=2)
#
#
#B = collect(transpose(atoms))
#
#@time CH = QHull.chull(B)
#
#
#CH = transpose(B[CH.vertices,:])
#
##	Plots.plot!(h[1,:],h[2,:])
#	Plots.plot!(CH[1,:],CH[2,:])
#
#
#
#
#	plot!(latts[end],xlims=xlim,ylims=ylim,ratio=1)
#
#
##	extrema.(PosAtoms.(latts),dims=2)
#
#
#
#
##	return map(latts) do latt
##
##		Rs = PosAtoms(latt)
##
##		return Plots.scatter!(Rs[1,:],Rs[2,:])
##
##	end 
##
#
#
#end 
#
#
#




#############################################################################
end
