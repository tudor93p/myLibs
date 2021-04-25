module Lattices
#############################################################################

using OrderedCollections:OrderedDict

import LinearAlgebra; const LA=LinearAlgebra
import Utils, Algebra

import Plots, ConcaveHull, QHull




export Lattice



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

const TOLERANCE = 1e-5


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

	LattVec::AbstractMatrix{Float64}
	
	Sublattices::OrderedDict{Any, AbstractMatrix{Float64}}

	Vacancies::OrderedDict{Any, AbstractMatrix{Float64}}

	function Lattice(LV, SL=nothing, VA=nothing; mode::Symbol=:cartesian)
									 
		lv, sl, va = to_myMatrix(LV), to_myODict(SL), to_myODict(VA) 
		
		for k in [keys(sl);keys(va)], odict in [sl,va]

			if haskey(odict, k)

				s = size(odict[k],1)

				if mode==:fractional

					size(lv,2)!=s && error("Dimensionality mismatch. $s coefficient(s) was(were) provided for sublattice '$k', but ",size(lv,2)," are necessary")

					odict[k] = CombsOfVecs(lv, odict[k])
				
				elseif mode==:cartesian

					size(lv,1)!=s && error("Dimensionality mismatch. The positions for sublattice '$k' are $s-dimensional, but the lattice vectors are ",size(lv,1),"-dimensional")

				else 

					error("Mode $mode not defined")

				end 


			end 

		end 

		return new(lv,sl,va)

	end 

	function Lattice(latt::Lattice; act_on_vectors=(v::AbstractMatrix)->v,
									 								act_on_atoms=(d::AbstractDict)->d
																								)
		Lattice(
	
			collect( 

				if hasmethod(act_on_vectors, (AbstractMatrix,))

					act_on_vectors(copy(latt.LattVec))

				elseif hasmethod(act_on_vectors, (Lattice,))

					act_on_vectors(latt)

				else

					error("Improper definiton of 'act_on_vectors'")

				end

						),

		
			map([:Sublattices, :Vacancies]) do kind
	
				D = getproperty(latt, kind)

				for args in ((D,), (kind,), (), (kind, D),)

					!hasmethod(act_on_atoms, typeof.(args)) && continue 
			
					return act_on_atoms(map(copy, args)...)

				end 

	
				return map(collect(keys(D))) do k 

					for args in ((D[k], k), (kind, D[k], k), (kind, D[k]))

						!hasmethod(act_on_atoms, typeof.(args)) && continue 

						return act_on_atoms(map(copy, args)...)

					end 

					error("Improper methods of 'act_on_atoms'")

				end 

			end...)
	
	
	end 
	
		
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

	error("Wrong dimensions or types")

end 


function to_myMatrix(x, latt::Lattice)::AbstractMatrix{Float64}

	to_myMatrix(x, VecDim(latt))
	
end 

function to_myMatrix(x, A::AbstractMatrix)::AbstractMatrix{Float64}

	to_myMatrix(x, size(A,1))
	
end 

#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#



function to_myODict(x::Nothing, D::Union{Int64,Nothing}=nothing
									 )::OrderedDict{Any, AbstractMatrix{Float64}}

	OrderedDict()

end


function to_myODict(x::AbstractDict, D::Union{Int64,Nothing}=nothing
									 )::OrderedDict{Any, AbstractMatrix{Float64}}

	if x isa OrderedDict 

		valtype(x)==AbstractMatrix{Float64} && return x

		return OrderedDict(k=>to_myMatrix(v, D) for (k,v) in x)

	end 

	return OrderedDict(map(sort(collect(keys(x)), by=string)) do k

									 k => to_myMatrix(x[k], D)

								 end)

end 


function to_myODict(x::Utils.List, D::Union{Int64,Nothing}=nothing
									 )::OrderedDict{Any, AbstractMatrix{Float64}}

	isempty(x) && return OrderedDict()

	all(isa.(x,Pair)) || return OrderedDict("A"=>to_myMatrix(x, D))
	
	return OrderedDict(k=>to_myMatrix(v, D) for (k,v) in x)

end


function to_myODict(x::AbstractMatrix, D::Union{Int64,Nothing}=nothing
									 )::OrderedDict{Any, AbstractMatrix{Float64}}

	OrderedDict("A"=>to_myMatrix(x, D))

end



function to_myODict(x::Pair, D::Union{Int64,Nothing}=nothing
									 )::OrderedDict{Any, AbstractMatrix{Float64}}

	OrderedDict(x.first => to_myMatrix(x.second, D))

end


function to_myODict(x, latt::Lattice
									 )::OrderedDict{Any, AbstractMatrix{Float64}}

	to_myODict(x, VecDim(latt))

end 

function to_myODict(x, A::AbstractMatrix
									 )::OrderedDict{Any, AbstractMatrix{Float64}}

	to_myODict(x, size(A,1))

end 

#to_myODict(x...)::OrderedDict{Any, AbstractMatrix{Float64}} = to_myODict(x)




#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function CombsOfVecs(A::AbstractMatrix{<:T}, coeff::AbstractMatrix{<:V}; dim=2)::AbstractMatrix{promote_type(T,V)} where T where V

	dim==2 && return A*coeff

	dim==1 && return coeff*A

	error()

end



function CombsOfVecs(latt::Lattice, coeff; kwargs...)::AbstractMatrix{Float64}

	CombsOfVecs(LattVec(latt), to_myMatrix(coeff, LattDim(latt)); kwargs...)

end






function CombsOfVecs(t::Tuple{<:Union{Lattice, AbstractMatrix{<:Real}},Any};
										 kwargs...)::AbstractMatrix{Float64}

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

function sublatt_labels(latt::Lattice; kind=:Sublattices)

	kind!=:Both && return collect(keys(getproperty(latt, kind)))

	return union(sublatt_labels(latt; kind=:Sublattices),
							 sublatt_labels(latt; kind=:Vacancies))

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


function sublattices_contain(latt::Lattice, inp::AbstractString=""; kwargs...)

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


EmptyPos(D::Int64)::AbstractMatrix{Float64} = zeros(Float64, D, 0)

EmptyPos(n::Nothing)::AbstractMatrix{Float64} = EmptyPos(0)

EmptyPos(latt::Lattice)::AbstractMatrix{Float64} = EmptyPos(VecDim(latt))


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#




function PosAtoms(latt::Lattice;
									labels_contain=nothing,
									label=nothing,
									f=(sl,P)->P,
									kind=:Sublattces,
									kwargs...
									)

	SL = if !isnothing(label)
		
					in(label, sublatt_labels(latt, kind=kind)) ? [label] : []
					
					else 
		
						sublattices_contain(latt, labels_contain)

				end 




	for sl in SL 

		isempty(latt.Sublattices[sl]) && continue

		return hcat([f(s, latt.Sublattices[s]) for s in SL]...)

	end 


	return EmptyPos(latt)
	

end


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function VecDim(latt::Lattice)::Int

	size(latt.LattVec,1)>0 && return size(latt.LattVec,1) 

	for prop in [:Sublattices, :Vacancies], (k,v) in getproperty(latt, prop)
		
		!isnothing(v) && size(v,1)>0 && return size(v,1)
	
	end 

	return 0

end 


function LattDim(latt::Lattice)::Int

	size(latt.LattVec,2)

end 


function LattVec(latt::Lattice)::AbstractMatrix{Float64}

	latt.LattVec

end 

LattVec(latt::Nothing)::Nothing = nothing


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function CombsOfVecs10(A_or_latt::T, stopstart...; dim=2) where T<:Union{Lattice,AbstractMatrix{<:Real}}

	D = if T<:Lattice 

					LattDim(A_or_latt) 
					
			elseif T<:AbstractMatrix

					size(A_or_latt, dim)

			else 

				error("Type '$T' not supported")

			end 

	return CombsOfVecs(A_or_latt,
										 Utils.vectors_of_integers(D, stopstart...; dim=dim),
										 dim=dim)

end 



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function UnitCells(latt::Lattice, stopstart...)

	CombsOfVecs10(latt, stopstart...; dim=2)

end 

#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function OuterDist(args::Vararg{Any,2}; kwargs...)::AbstractArray

	Algebra.OuterDist(map(args) do X

		isa(X,Lattice) ? PosAtoms(X; kwargs...) : X

	end...; dim=2)

end 


function FlatOuterDist(args::Vararg{Any,2}; kwargs...)::AbstractArray

	Algebra.FlatOuterDist(map(args) do X

		isa(X,Lattice) ? PosAtoms(X; kwargs...) : X

	end...; dim=2)

end 







function Distances(latt::Lattice; nr_uc=2, nr_neighbors=1, with_zero=true, kwargs...)

	TOLERANCE>1e-1 && error("The tolerance is too large!")

	AtomsUC = PosAtoms(latt; kwargs...)
	
	Ds = FlatOuterDist(AtomsUC[:, 1:1].-AtomsUC, UnitCells(latt, nr_uc))

	Utils.Unique!(Ds; tol=TOLERANCE, sorted=true)

	!with_zero && deleteat!(Ds, Ds.<TOLERANCE)

	return nr_neighbors=="all" ? Ds : Ds[1:min(nr_neighbors+1,end)]
  
end 




#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#



function SquareInt_LattCoeff((latt,n)::Tuple{Lattice, T})::AbstractMatrix{Int} where T

	SquareInt_LattCoeff(latt, n)

end
														

function SquareInt_LattCoeff(latt::Lattice, n)::AbstractMatrix{Int}

	n = to_myMatrix(n, latt)

	tn = trunc.(n)

	if !isapprox(tn, n, atol=TOLERANCE) 

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
# Generate vertices of a d-dimensional body based on d vectors
#
#---------------------------------------------------------------------------#

function BodyVertices_fromVectors(v::AbstractMatrix{T}; dim=2)::AbstractMatrix{T} where T

	CombsOfVecs10(v, 1, 0; dim=dim)

end 



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function is_left(P0, P1, P2)

	(P1[1] - P0[1]) * (P2[2] - P0[2]) - (P2[1] - P0[1]) * (P1[2] - P0[2])

end 

function Order_PolygonVertices(V::AbstractMatrix; dim=2)

	C = sum(V, dims=dim)[:]/size(V,dim) # center of the polygon

	return sortslices(V, dims=dim, by = R->atan(R[2]-C[2], R[1]-C[1]))

end 


function prepare_polygon_vertices(V::AbstractMatrix; 
																	order_vertices=false, dim=2)::AbstractMatrix

	if order_vertices

		return prepare_polygon_vertices(
								Order_PolygonVertices(V, dim=dim),
								order_vertices=false, dim=dim)

	end 

	return cat(selectdim(V, [2,1][dim], 1:2),
					selectdim(selectdim(V, [2,1][dim], 1:2), dim, 1:1),
					dims=dim)

end 

function PointInPolygon_wn(V::AbstractMatrix; kwargs...)

	V = prepare_polygon_vertices(V; kwargs...) 
															 
	return P -> PointInPolygon_wn(P, V; kwargs..., prepare_vertices=false)

end 



function PointInPolygon_wn(P::Utils.List, V::AbstractMatrix; 
													 dim=2, prepare_vertices=true, kwargs...)

  wn = 0   # the winding number counter

  						# repeat the first vertex at the end

	if prepare_vertices

		V = prepare_polygon_vertices(V; dim=dim, kwargs...)

	end 

	v(i) = selectdim(V, dim, i)

	v(i,j) = v(i)[j]



  # loop through all edges of the polygon

	for i in 1:size(V, dim)-1  # edge from v(i) to v(i+1) 

		if v(i,2) <= P[2] # start y <= P[1]

			if v(i+1,2)>P[2] && is_left(v(i), v(i+1),P)>0 
									# an upward crossing & P left of edge
			
				wn += 1           # have a valid up intersect

			end 

			 
		else # start y > P[2] (no test needed)
			 
			if v(i+1,2)<=P[2] && is_left(v(i), v(i+1), P)<0

						
						# a downward crossing &  P right of edge

				wn -= 1           # have a valid down intersect


			end

		end 

	end 

  return wn != 0


end 




#===========================================================================#
#
# Obtain the ucs in big UC by checking which ucs are inside a polygon
#
#---------------------------------------------------------------------------#

function ucsUC_Polygon(v, mM, polygon; dim=2, asfunction=false)

	is_inside = PointInPolygon_wn(polygon, order_vertices=true, dim=dim)


	uc_candidates = Utils.vectors_of_integers(2, 	selectdim(mM, dim, 2).+1,
																								selectdim(mM, dim, 1).-1;
																								dim=dim)


	out(cand) = selectdim(cand, dim, 
												mapslices(is_inside, cand, dims=[2,1][dim])[:]
												) #|> collect


	!asfunction && return out(uc_candidates)


	shifts = sortslices(CombsOfVecs10(polygon, 1; dim=dim),
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

function ucs_in_UC(N::AbstractMatrix{T}; dim=2, method2D_cells=:Polygon, kwargs...) where T

	if T<:Int 

		# Nothing to do, N is ready.

	elseif T<:Float64 || any(n->isa(n,Float64), N)

		tN = trunc.(N) 

		isapprox(tN, N, atol=TOLERANCE) || error("N should contain integers")
		
		return ucs_in_UC(convert(AbstractMatrix{Int}, tN); kwargs...)

	else 

		error("Type '$T' not supported")

	end 


	correct_nr = Int(round(abs(LA.det(N))))
	# this should be the size of the new unit cell -- like | a1.(a2 x a3) |


	polygon = BodyVertices_fromVectors(N; dim=dim) 

	mM = cat(cat.(extrema(polygon, dims=dim)..., dims=[2,1][dim])..., dims=dim)


#  if len(N) == 1:
#
#    iter_ucs = [np.arange(*maxmin[::-1]).reshape(-1,1)]
#

	n, iter_ucs = if size(N,2)==2 && method2D_cells == :Polygon
									
									ucsUC_Polygon(N, mM, polygon; dim=dim, asfunction=true)

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

		size(ucs, dim)==correct_nr && return ucs
		
		# if the correct number of uc-s in the large UC is found

	end 

	error("Could not construct big unit cell correctly.")

end

#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

#function CollectAtoms(lattices::Utils.List, latt_vac::Symbol, label)::AbstractMatrix{Float64}
#
#	hcat(map(lattices) do L 
#
#		L::Lattice
#
#		!hasproperty(L, latt_vac) && error("'$latt_vac' not found")
#
#		return get(getproperty(L, latt_vac), k, zeros(Float64, VecDim(L), 0))
#
#	end...) 
#
#end 

function parse_input_RsNs( latt::Union{Nothing,Lattice}=nothing;
													 Rs=nothing,
													 Ns=nothing,
													 A=nothing,
													 StopStart=nothing,
													 dim=2,
													 kwargs...
													)

	!isnothing(Rs) && return Rs 

	A = Utils.Assign_Value(A, LattVec(latt))

	isnothing(A) && error("Cannot proceed with unknown vectors")

	!isnothing(Ns) && return CombsOfVecs(A, to_myMatrix(Ns, latt), dim=dim)

	!isnothing(StopStart) && return CombsOfVecs10(A, StopStart...; dim=dim)

	error("Not enough non-nothing kwargs")

end




#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#




function Atoms_ManyUCs(atoms::AbstractMatrix{<:Real},
											 ucs::AbstractMatrix{<:Real};
											 dim::Int=2, kwargs...)::AbstractMatrix{Float64}

	Algebra.FlatOuterSum(atoms, ucs; dim=dim)

end 


function Atoms_ManyUCs(atoms::AbstractMatrix{<:Real},
											 ns::AbstractMatrix{<:Real},
											 A::AbstractMatrix{<:Real};
											 kwargs...)::AbstractMatrix{Float64}

	Atoms_ManyUCs(atoms, parse_input_RsNs(;Ns=ns, A=A, kwargs...))

end 



function Atoms_ManyUCs(latt::Lattice; kwargs...)::AbstractMatrix{Float64}

	Atoms_ManyUCs(PosAtoms(latt; kwargs...), parse_input_RsNs(latt; kwargs...))

end







#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function Superlattice(latt::Lattice, n; kwargs...)
	
	Superlattice([latt], [n]; kwargs...)

end





function Superlattice(Components::Utils.List, Ns::Utils.List; Labels=nothing, kwargs...)

	Ns = SquareInt_LattCoeff.(zip(Components,Ns))

	Supervectors = CombsOfVecs(Components[1], Ns[1])

	for S in CombsOfVecs.(zip(Components[2:end], Ns[2:end]))

		isapprox(Supervectors, S, atol=TOLERANCE) && continue

		error("The 'Ns' provided are not compatible; one should recover the same supervectors for each pair (N,Component).")
	
	end



	Ns = ucs_in_UC.(Ns; kwargs...)

	Labels = 	if isnothing(Labels)
	
							repeat([""], length(Components))

						else 

							Utils.Assign_Value.(Labels, "")

						end 

	empty_labels = isempty.(Labels)


	return Lattice(Supervectors, map([:Sublattices, :Vacancies]) do p

		Utils.flatmap(unique(vcat(sublatt_labels.(Components; kind=p)...))) do k 

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

	end...) # map over [:Sublattices, :Vacancies] ends

end





function Superlattice!(latt::Lattice, n; kwargs...)

	N = SquareInt_LattCoeff(latt, n)
	
	ucs = ucs_in_UC(N; kwargs...)

	for p in [:Sublattices, :Vacancies], k in sublatt_labels(latt; kind=p)

		getproperty(latt, p)[k] = Atoms_ManyUCs(latt; Ns=ucs, label=k, kind=p)

	end

	latt.LattVec = CombsOfVecs(latt, N)

	return latt

end



#===========================================================================#
#
# 			some safety measureas when adding atoms 
#
#
#---------------------------------------------------------------------------#

function Initialize_PosAtoms(latt::Lattice, pos_atoms...)

	Initialize_PosAtoms(VecDim(latt), pos_atoms...)

end 

function Initialize_PosAtoms(vect_dim::Int, pos_atoms=nothing)
	
#    if type(pos_atoms)==type(None) or np.size(pos_atoms)==0:
#      print("\n *** No position given. The atom is set at the origin. ***\n")

	# if nothing is given, take the origin

	pos_atoms = to_myMatrix(if isnothing(pos_atoms) || isempty(pos_atoms) 

														zeros(Float64, vect_dim, 1)

													else 

														pos_atoms

													end)

	dim,nr = size(pos_atoms)

	if dim>vect_dim

		error("Atoms cannot have more coordinates than the dimension of the lattice vectors!")

	elseif dim==vect_dim

		return pos_atoms

	else 

	# in case additional coordinates are needed, use zeross

		return vcat(pos_atoms, zeros(Float64, vect_dim-dim, nr))

	end 

end 



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function parse_input_sublattices(latt::Lattice, inp::Nothing, args...)

	good = !in(sublatt_labels(latt; kind=:Both))

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



function AddAtoms!(latt::Lattice, positions=[], labels="A"; kind=:Sublattices)

	positions = Initialize_PosAtoms(latt, positions)

	new_labels = parse_input_sublattices(latt, labels, size(positions,2), string)

	D = getproperty(latt, kind)


	for (k,i) in zip(Utils.Unique(new_labels, inds=:all)...)

		isempty(i) && continue

		D[k] = hcat(get(D, k, EmptyPos(latt)), view(positions, :, i))

	end 

	return latt 

end 


function AddAtoms(latt::Lattice, positions=[], labels="A"; kind=:Sublattices)

	positions = Initialize_PosAtoms(latt, positions)

	new_labels, inds = Utils.Unique(
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

function ShiftAtoms!(latt::Lattice; n=nothing, r=nothing, kind=:Both)

	R = parse_input_RsNs(latt; Ns=n, Rs=r, dim=2)

	for K in (kind==:Both ? [:Sublattices, :Vacancies] : [kind])

		for l in sublatt_labels(latt, kind=K)
			
			getproperty(latt, K)[l] .+= R 

		end 

	end 

	return latt 

end


function ShiftAtoms(latt::Lattice; n=nothing, r=nothing, kind=:Both)

	R = parse_input_RsNs(latt; Ns=n, Rs=r, dim=2)

	shift(v::AbstractMatrix{Float64}, l) = Pair(l, v.+R)


	kind==:Both && return Lattice(latt, act_on_atoms=shift)


	function act_on_atoms(K::Symbol, D::AbstractDict)
	
		kind==K ? [shift(v,l) for (l,v) in D] : D

	end 

	return Lattice(latt, act_on_atoms=act_on_atoms)


end 




#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function ReduceDim(latt::Lattice, dim=nothing; complement=false)

	if complement && (isnothing(dim) || length(dim)>1)
		
		error("The complement can be computed only if one single dimension is given.")

	end 

	complement && return Lattice(latt, act_on_vectors=v->v[:, dim[1]:dim[1]])

	isnothing(dim) && return Lattice(latt, act_on_vectors=EmptyPos)

	keep_dims = filter(!in(vcat(dim...)),1:LattDim(latt))

	return Lattice(latt, act_on_vectors=v->v[:, keep_dims])

end 


   
function ReduceDim!(latt::Lattice, dim=nothing; complement=false)

	if complement && (isnothing(dim) || length(dim)>1)
		
		error("The complement can be computed only if one single dimension is given.")

	end 


	keep_dims = if complement 

					[dim[1]]

							else 

					isnothing(dim) ? [] : filter(!in(vcat(dim...)),1:LattDim(latt))

						end

	latt.LattVec = latt.LattVec[:, keep_dims]

	return latt 

end 




#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function SurfaceAtoms(latt::Lattice, bondlength=nothing; kwargs...)

	LattDim(Latt)==0 || error()

	SurfaceAtoms(PosAtoms(latt; kwargs...), bondlength)

end 


function SurfaceAtoms(atoms::AbstractMatrix{Float64})

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


end 



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

function Align_toAtoms(latt::Lattice, Atoms::AbstractMatrix{Float64}, shift_dir::Int=+1; bonds=nothing)

	# shift_dir = +/- i means positive/negative a_i
	# slowed down if all atoms provided instead of the surface ones

	if isnothing(bonds)
		
		c = Utils.DistributeBallsToBoxes.([-1,1], LattDim(latt))

		bonds = CombsOfVecs(latt, vcat(c...))

	end 

	Atoms = Order_PolygonVertices(Atoms)


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



	dmax = maximum(Algebra.OuterDist(Atoms, Atoms))*5


	dmax=10

	direction = sign(shift_dir)*LattVec(latt)[:,abs(shift_dir)]


	longdir = 2*dmax * LA.normalize(direction)


#	ShiftAtoms!(latt, r=-longdir/2)


	Maximize_ContactSurface(PosAtoms(latt),
													longdir,
													ordered_vertices=Atoms)

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
import 2libs 

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








#############################################################################
end
