module Operators
#############################################################################




import ..LA, ..SpA 

import ..Utils, ..TBmodel, ..Algebra



#===========================================================================#
#
# Basic sturcture
#
#---------------------------------------------------------------------------#



struct Operator 

	data::Union{Number, 
							<:VecOrMat{<:Number},
							<:Vector{<:VecOrMat{<:Number}}}

	diag::Symbol  


#	nr_at::Int # remove, eventually
#	
#	nr_orb::Int # remove 
#	
#	size_H::Int # remove 


	vsdim::Int 

	csdim::Int

#	selvect::Function 
#
#	selcomp::Function



	sums::NTuple{2,Bool}


#	acts_on_atoms::Vector{Int}

#	acts_on_orbs::Vector{Int}


	inds::Union{Vector{Int}, Vector{Vector{Int}}, Colon}

#	iter_Op::Any





	function Operator(; kwargs...)

		Operator(1; kwargs...)

	end

	function Operator(diag::Symbol; kwargs...)

		Operator(1, diag; kwargs...)

	end


	function Operator(A::Union{Number, AbstractVecOrMat{<:Number}}, 
										args...; kwargs...)

		numbers,enough_data = get_nratorbH(; kwargs...)

		@assert enough_data "Provide more size-related kwargs"
	
		return Operator(reduce_data(A), numbers, args...; kwargs...)

	end   


	function Operator(data_::Union{Number, AbstractVecOrMat{<:Number}},
										numbers::NTuple{3,Int},
										diag...;
										kwargs...)

		diag = prep_diag(data_, numbers, diag...)

		ws = which_sum(;kwargs...)

		ao = get_acts_on(diag, numbers...; kwargs...)
	
		dim = get_vc_dims(;kwargs...)


		data = extend_data(data_, Val(diag), get_nratorbH(length.(ao)...)[1], Val.(ws))

		inds = get_iter_inds(ws, ao, numbers) 

		iter_Op = get_iter_Op(data, Val.(ws), inds, Val(diag))

		return new(iter_Op, diag, dim..., ws, inds)

	end  

end 




#===========================================================================#
#
# 
# 
#---------------------------------------------------------------------------#

function reduce_data(A::AbstractMatrix{T})::Union{T,VecOrMat{T}} where T<:Number 

	size(A,1) != size(A,2) && return reduce_data(reshape(A,:))

	dA = LA.diag(A)

	isapprox(A, LA.diagm(0=>dA), atol=1e-6) && return reduce_data(dA)

	return A 

end  

function reduce_data(A::AbstractVector{T})::Union{T,Vector{T}} where T<:Number 

	length(Utils.Unique(A; tol=1e-6))==1 && return A[1]
	
	return Vector{T}(A)

end 

function reduce_data(A::T)::T where T<:Number 

	A

end


function prep_diag(A::Union{Number,<:AbstractVecOrMat{<:Number}},
									 numbers::NTuple{3,Int})::Symbol

	possib_diag = which_diag(A, numbers)

	if :atoms in possib_diag && :orbitals in possib_diag 

		@assert numbers[1]==numbers[2]==1 "Confusing sizes. Provide 'diag' arg."

	end  

	return first(possib_diag)

end 


function prep_diag(A::Union{Number,<:AbstractVecOrMat{<:Number}},
									 numbers::NTuple{3,Int}, diag::Symbol)::Symbol

	possib_diag = which_diag(A, numbers) 

	for (weak,strong) in [([diag],[diag]),
												([:none],[:all,:orbitals,:atoms]),
												([:orbitals,:atoms],[:all])]

		diag in weak || continue 

		for k in strong 

			k in possib_diag && return k

		end 

	end 



	error("The 'diag' given is not consistent with sizes")


end 
















function get_nratorbH(nr_at::Integer, nr_orb::Integer, size_H::Integer
											)::Tuple{NTuple{3,Int},Bool}

	(nr_at, nr_orb, size_H), nr_at*nr_orb==size_H

end


function get_nratorbH(nr_at::Nothing, nr_orb::Integer, size_H::Integer
											)::Tuple{NTuple{3,Int},Bool}

	get_nratorbH(div(size_H,nr_orb), nr_orb, size_H)

end

function get_nratorbH(nr_at::Integer, nr_orb::Nothing, size_H::Integer
											)::Tuple{NTuple{3,Int},Bool}

	get_nratorbH(nr_at, div(size_H,nr_at), size_H)

end

function get_nratorbH(nr_at::Integer, nr_orb::Integer, size_H::Nothing=nothing
											)::Tuple{NTuple{3,Int},Bool}

	get_nratorbH(nr_at, nr_orb, nr_at*nr_orb)

end


function get_nratorbH(;nr_at=nothing, nr_orb=nothing, size_H=nothing, 
											kwargs...)::Tuple{NTuple{3,Int},Bool}

	args = (nr_at,nr_orb,size_H) 

	count(a->isa(a,Int),args)>=2 && return get_nratorbH(args...)

	return Tuple(isa(a,Int) ? a : 0 for a in args), false 

end




#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


which_diag(A::Number, args...)::Vector{Symbol} = which_diag([A], args...)

function which_diag(A::AbstractVecOrMat, nrs::NTuple{3,Int})::Vector{Symbol}

	out = Symbol[]

	prod(size(A))==1 && push!(out, :all) 
	
	# these are worse (more expensive) solutions and come after ':all'
	# the better sols come in front  


	for k in [:orbitals, :atoms, :none][reverse(findall(nrs.==size(A,1)))]

		push!(out, k) 

	end  


	@assert !isempty(out) "Wrong operators size" 

	return out 

end 


function which_sum(trace::Vector{Symbol})::NTuple{2,Bool}

	K = [:atoms, :orbitals]

	@assert isempty(setdiff(trace, K))

	return Tuple(in.(K,[trace]))

end 

function which_sum(trace::Symbol)::NTuple{2,Bool} 
	
	trace==:none && return (false, false)
	trace==:all && return (true, true)
	
	return which_sum([trace])

end



function which_sum(;dim::Int, kwargs...)::Tuple{Bool,Bool}

	K = [:sum_atoms, :sum_orbitals]

	atorb = haskey.([kwargs], K) 

	D = (dim,[2,1][dim])

	if haskey(kwargs, :trace) 
	
		any(atorb) && error()

		return which_sum(kwargs[:trace])

	end 

	return Tuple((hk ? kwargs[k] : true) for (hk,k) in zip(atorb,K))

end

function get_vc_dims(;dim::Int, kwargs...)::Tuple{Int,Int}#,Function,Function}

	dim==1 && @warn "'dim=1' will most likely not work well!"

	(dim, [2,1][dim])#, Utils.sel(dim), Utils.sel([2,1][dim]))

end 


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function get_acts_on(diag::Symbol, 
										 nr_at::Int, nr_orb::Int, size_H::Int;
										 kwargs...)::NTuple{2,Vector{Int}}

	(vcat(get(kwargs, :acts_on_atoms, 1:nr_at)...),
	 vcat(get(kwargs, :acts_on_orbs, 1:nr_orb)...))


end

#===========================================================================#
#
# Compute expectation value
#
#---------------------------------------------------------------------------#


function get_iter_inds(sums::NTuple{2,Bool}, args...
											 )::Union{Colon, Vector{Int}, Vector{Vector{Int}}}

	get_iter_inds(Val.(sums), args...)

end 


function get_iter_inds(::Tuple{Val{true}, Val{false}},
											 (atoms,orbs)::NTuple{2,<:AbstractVector{Int}},
											 (nr_at,nr_orb,size_H)::NTuple{3,Int}
											 )::Vector{Vector{Int}}

	TBmodel.Hamilt_indices_all(orbs, atoms, nr_orb; iter="orbitals")

end 


function get_iter_inds(::Tuple{Val{false}, Val{true}},
											 (atoms, orbs)::NTuple{2,<:AbstractVector{Int}},
											 (nr_at,nr_orb,size_H)::NTuple{3,Int}
											 )::Vector{Vector{Int}}

	TBmodel.Hamilt_indices_all(orbs, atoms, nr_orb; iter="atoms") 

end  


function get_iter_inds(::NTuple{2,T},
											 acts_on::NTuple{2,<:AbstractVector{Int}},
											 numbers::NTuple{3,Int}
											 )::Union{Colon, Vector{Int}
																} where T<:Union{Val{true}, Val{false}}

	L = length.(acts_on)

	all(L.==length.(numbers[1:2])) && return Colon()

	return sort!(vcat(get_iter_inds(argmax(L).==(1,2), acts_on, numbers)...))

end 


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function check_size(A::AbstractVecOrMat{<:Number},
										sums::NTuple{2,Val},
										I::AbstractVector{<:AbstractVector{Int}},
										args...)

	check_size(A, sums,
						 get_nratorbH(length(I), only(unique(length.(I))))[1],
						 args...)

end  


function check_size(A::AbstractVecOrMat{<:Number},
										::Tuple{Val{true}, Val{false}},
										(L1,L2,L12)::NTuple{3,Int},
										diag::Val{:orbitals}
										)

#	println("\n***\n")
#	@show size(A) L1 L2 diag 

	@assert only(unique(size(A)))==L2

end 

function check_size(A::AbstractVector{<:Number},
										::Tuple{Val{true}, Val{false}},
									 (L1,L2,L12)::NTuple{3,Int},
										::Val{:atoms}
										)

	@assert only(unique(size(A)))==L1

end 


function check_size(A::AbstractMatrix{<:Number},
										sums::Tuple{Val{true}, Val{false}},
										numbers::NTuple{3,Int},
										diag::Val{:atoms}
										)
	@show size(A) sums diag 

	error("Op::Matrix => only '(false,true)+atoms' and '(true,false)+orbitals are allowed.")

end 

function check_size(A::AbstractVecOrMat{<:Number},
										sums::Tuple{Val{false},Val{true}},
										numbers::NTuple{3,Int},
										diag::Union{Val{:atoms},Val{:orbitals}})

	check_size(A, reverse(sums), numbers,
						 only(setdiff([Val(:atoms),Val(:orbitals)],[diag]))
						 )
end 


function check_size(A::AbstractVecOrMat{<:Number},
										::NTuple{2,Val},
										(L1,L2,L12)::NTuple{3,Int},
										::Val{:none}
										)

	@assert only(unique(size(A)))==L12 

end 

function check_size(A::Number, ::NTuple{2,Val}, arg, ::Val{:all})

end


function check_size(A::AbstractVecOrMat{<:Number},
										::NTuple{2,S} where S<:Union{Val{true},Val{false}},
										I::AbstractVector{Int},
										args...)

	@assert only(unique(size(A)))==length(I)

end 

#function check_size(A::AbstractVecOrMat,
#										sums::NTuple{2,S} where S<:Union{Val{true},Val{false}},
#										 (L1,L2,L12)::NTuple{3,Int},
#										 diag::Union{Val{:orbitals},Val{:atoms}}
#										) 
#	println("\n***\n")
#
#	@show size(A) sums L1 L2 L12  diag 
#
#	@assert only(unique(size(A)))==L12 

#end 


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function get_iter_Op(args...)
	
	check_size(args...)

	return get_iter_Op_(args...)

end 



function get_iter_Op_(Op::T,
									 sums::Tuple{Val{true}, Val{false}},
									 I::AbstractVector{<:AbstractVector{Int}},
									 diag::Val{:atoms}
									 )::T where T<:AbstractVecOrMat{<:Number}

	Op

end 


function get_iter_Op_(Op::T,
									 sums::Tuple{Val{true}, Val{false}},
									 I::AbstractVector{<:AbstractVector{Int}},
									 diag::Val{:orbitals},
									 )::Vector{T} where T<:AbstractVecOrMat{<:Number}

	[Op for i in I]

end 




function get_iter_Op_(Op::T,
									 sums::Union{Tuple{Val{true}, Val{false}},
													 Tuple{Val{false}, Val{true}}},
									 I::AbstractVector{<:AbstractVector{Int}},
									 diag::Val{:none},
									 )::Vector{T} where T<:AbstractVecOrMat{<:Number}

	[Op[(i for _=1:ndims(Op))...] for i in I]

end  


function get_iter_Op_(Op::T,
										 sums::Union{Tuple{Val{true},Val{false}},
														 Tuple{Val{false},Val{true}}},
									 I::AbstractVector{<:AbstractVector{Int}},
									 diag::Val{:all})::Vector{T} where T<:Number 

	[Op for i in I]

end 




function get_iter_Op(Op::T,
									 sums::Tuple{Val{false}, Val{true}},
									 I::AbstractVector{<:AbstractVector{Int}},
									 diag::Union{Val{:atoms},Val{:orbitals}},
									 )::Union{T,Vector{T}} where T<:AbstractVecOrMat{<:Number}

	get_iter_Op(Op, reverse(sums), I,
							only(setdiff([Val(:atoms),Val(:orbitals)],[diag])))

end 


function get_iter_Op_(Op::T,
											sums::NTuple{2,S},
											args...)::T where {T,S<:Union{Val{true},Val{false}}}
	Op 

end 


#===========================================================================#
#
# select only a subset of the wavefunctions
#
#---------------------------------------------------------------------------#


function ExpectVal(P::AbstractMatrix{<:Number}, 
									 Op::Number,
									 sums::NTuple{2,Val{false}},
									 csdim::Int,
									 I::AbstractVector{Int},
									 args...)::Matrix{<:Number}

	ExpectVal(selectdim(P, csdim, I), Op, sums, csdim, :, args...) 

end 




function ExpectVal(P::AbstractMatrix{<:Number}, 
									 Op::AbstractMatrix{<:Number},
									 sums::NTuple{2,Val{true}},
									 csdim::Int,
									 I::AbstractVector{Int},
									 args...)::Matrix{<:Number}

	ExpectVal(selectdim(P, csdim, I), Op, sums, csdim, :, args...)

end 




#===========================================================================#
#
# Basic methods for computing the expectation values 
# 								(in the simple cases sum all or nothing)
#
#---------------------------------------------------------------------------#

function ExpectVal(P::AbstractMatrix{<:Number}, 
									 Op::Number,
									 ::NTuple{2,Val{false}},
									 csdim::Int,
									 I::Colon,
									 args...)::Matrix{<:Number}
	abs2.(P)*Op 

end 


function ExpectVal(P::AbstractMatrix{<:Number}, 
									 Op::AbstractVector{<:Number},
									 sums::NTuple{2,Val{false}},
									 csdim::Int,
									 args...)::Matrix{<:Number}

	aP = ExpectVal(P, 1, sums, csdim, args...)

	csdim==1 && return Op .* aP 
	
	csdim==2 && return aP .* Op 

end 


function ExpectVal(P::AbstractMatrix{<:Number}, 
									 Op::Union{Number,<:AbstractVector{<:Number}},
									 sums::NTuple{2,Val{true}},
									 csdim::Int,
									 args...)::Matrix{<:Number}

	sum(ExpectVal(P, Op, (Val(false),Val(false)), csdim, args...), dims=csdim)

end  


function ExpectVal(P::AbstractMatrix{<:Number}, 
									 Op::AbstractMatrix{<:Number},
									 ::NTuple{2,Val{true}},
									 csdim::Int,
									 I::Colon,
									 args...)::Matrix

	Utils.VecAsMat(LA.diag(P' * Op * P), csdim)

end 




#===========================================================================#
#
# Build matrix when only one of atoms/orbitals is iterated 
#
#---------------------------------------------------------------------------#


function ExpectVal(P::AbstractMatrix{<:Number},
										iter_Op::Union{T, <:AbstractVector{T}},
										::Union{Tuple{Val{true},Val{false}},
														Tuple{Val{false},Val{true}}},
									 csdim::Int,
									 I::AbstractVector{<:AbstractVector{Int}},
									 args...
									 )::Matrix{<:Number} where T<:AbstractVecOrMat{<:Number}

	Utils.VecsToMat(iter_Op, I; dim=csdim) do Op,inds

		ExpectVal(P, Op, (Val(true),Val(true)), csdim, inds)

	end 

end   

function Trace(A::AbstractVector{<:Number},
							 iter_Op::Union{T,<:AbstractVector{T}},
							I::AbstractVector{<:AbstractVector{Int}},
									 )::Vector{<:Number} where T<:AbstractVector{<:Number}

	[Trace(A[inds], Op) for (Op,inds) in zip(iter_Op, I)]

end 



function Trace(A::AbstractVector{<:Number},
							 Op::Union{Number, <:AbstractVector{<:Number}}=1,
							 I::Colon=Colon(),
							 )::Number 

	sum(A.*Op)

end 

function Trace(A::AbstractVector{<:Number},
							Op::Union{Number, <:AbstractVector{<:Number}},
							I::AbstractVector{Int},
							)::Number 

	Trace((A.*Op)[I])

end 


#===========================================================================#
#
# using an Operator as a one-argument function => compute expectation 
#
#---------------------------------------------------------------------------#

function (H::Operator)(P::AbstractMatrix{<:Number}; kwargs...
											)::Matrix{<:Number}

	ExpectVal(P, H.data, Val.(H.sums), H.csdim, H.inds, Val(H.diag))

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


extend_data(A::Number, args...)::Number = A



#function extend_data(A::T,
#										 diag::Val,
#										 numbers::NTuple{3,Int},
#										 sums::NTuple{2,Val}
#										 )::T where T<:AbstractVecOrMat{<:Number}
#
#
#	extend_data_(A, diag, numbers, sums)
#
##	check_size(data, sums, numbers, diag) 
#
##	return data 
#
#end



function extend_data(A::T, ::Val{:none}, args...
										 )::T where T<:AbstractVecOrMat{<:Number}
	A

end

function extend_data(A::T,
										 ::Union{Val{:atoms},Val{:orbitals}},
										 ::NTuple{3,Int},
										 ::Union{Tuple{Val{true},Val{false}},
															 Tuple{Val{false},Val{true}}},
										 )::T where T<:AbstractVecOrMat{<:Number}
	A

end 


function extend_data(A::T,
										 ::Val{:atoms},
										 numbers::NTuple{3,Int},
										 ::NTuple{2,V} where V<:Union{Val{true}, Val{false}}
										 )::T where T<:AbstractVecOrMat{<:Number}

	repeatOperator_manyAtoms(A, numbers...)

end




function extend_data(A::T,
										 ::Val{:orbitals},
										 numbers::NTuple{3,Int},
										 ::NTuple{2,V} where V<:Union{Val{true},Val{false}},
										 )::T where T<:AbstractVecOrMat{<:Number} 

	repeatOperator_manyOrbitals(A, numbers...)

end



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#



function repeatOperator_manyAtoms(Op::AbstractVector{T}, 
																	nr_at::Int, nr_orb::Int, size_H::Int
																	)::Vector{T} where T<:Number

	# Op simulates a Hamiltonian with nr_orb orbitals and one atom
	# 	will be reproduced for nr_at atoms

	@assert length(Op)==nr_orb 

	FullOp = similar(Op, size_H)

	for a in 1:nr_at 

		FullOp[TBmodel.Hamilt_indices(1:nr_orb, a, nr_orb)] = Op 

	end 

	return FullOp

end


function repeatOperator_manyAtoms(Op::AbstractMatrix{T},
																	nr_at::Int, nr_orb::Int, size_H::Int
																	)::Matrix{T} where T<:Number

	# Op simulates a Hamiltonian with nr_orb orbitals and one atom
	# 	will be reproduced for nr_at atoms

	@assert all(isequal(nr_orb), size(Op))
	
	FullOp = zeros(T, size_H, size_H)

	for a in 1:nr_at 

		i = TBmodel.Hamilt_indices(1:nr_orb, a, nr_orb) 

		FullOp[i,i] = Op 

	end 

	return FullOp

end






#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#




function repeatOperator_manyOrbitals(Op::AbstractVector{T}, 
																	nr_at::Int, nr_orb::Int, size_H::Int
																	)::Vector{T} where T<:Number

	# Op simulates a Hamiltonian with nr_at atoms and one orbital
	# 	will be reproduced for nr_orb orbitals

	@assert length(Op)==nr_at 

	FullOp = similar(Op, size_H)

	for orb in 1:nr_orb 

		FullOp[TBmodel.Hamilt_indices(orb, 1:nr_at, nr_orb)] = Op 

	end 

	return FullOp

end


function repeatOperator_manyOrbitals(Op::AbstractMatrix{T},
																	nr_at::Int, nr_orb::Int, size_H::Int
																	)::Matrix{T} where T<:Number

	@assert all(isequal(nr_at), size(Op))
	
	FullOp = zeros(T, size_H, size_H)

	for orb in 1:nr_orb

		i = TBmodel.Hamilt_indices(orb, 1:nr_at, nr_orb)

		FullOp[i,i] = Op 

	end 

	return FullOp

end



#===========================================================================#
#
# Trace over orbitals/atoms with some projection operator
# 
#---------------------------------------------------------------------------#


function Trace(what::AbstractString, args...; kwargs...)::Function 

	Trace(Symbol(what), args...; kwargs...)

end 

function Trace(what::Symbol,
								data::Union{Number,<:AbstractVecOrMat{<:Number}}=1;
								dim::Int,
								sum_up::Bool=false, kwargs...)::Function

	@assert what in [:orbitals, :atoms]

	function get_diag(X::AbstractMatrix{T})::Vector{T} where T<:Number

		LA.diag(X)

	end 

	function get_diag(X::T)::T where T<:Union{Number,AbstractVector{<:Number}}

		X  

	end 


	enough_data = get_nratorbH(; kwargs...)[2]

	data_ = get_diag(data)



	if !enough_data 
		
		@assert !haskey(kwargs, :size_H) 

		@assert haskey(kwargs, :nr_at) || haskey(kwargs, :nr_orb)

		return function tr2(A::AbstractVecOrMat{T}
												)::Union{T,Vector{T}} where T<:Number 

			a = get_diag(A)
		
			Op = Operator(data_; dim=dim, kwargs..., trace=sum_up ? :all : what,
										size_H=length(a))

			return Trace(a, Op.data, Op.inds)

		end 


	else 

		Op = Operator(data_; dim=dim, kwargs..., trace=sum_up ? :all : what)  


		return function tr(A::AbstractVecOrMat{T}
											 )::Union{T,Vector{T}} where T<:Number 
	
			Trace(get_diag(A), Op.data, Op.inds)

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
# Special operators usually used
# 
#---------------------------------------------------------------------------#



# ----- probability of localization on individual atoms ----- #


function LDOS(Op...; kwargs...)::Operator
	
	Operator(Op...; trace=:orbitals, kwargs...)

end 

function Position(ax::Int, Rs::AbstractMatrix{<:Real}; 
									fpos::Function=identity,
									dim::Int, nr_at::Int=size(Rs,dim),
									kwargs...)::Operator

	a = selectdim(Rs, [2,1][dim], ax)

	Op = applicable(fpos, a) ? fpos(a) : fpos.(a)

	return Operator(Op, :orbitals; dim=dim, nr_at=nr_at, kwargs...)

end 







	# ----- local charge operator ----- #


function SingleAtomCharge(Op_UC::AbstractVecOrMat{<:Number}, atom::Int; 
								kwargs...)::Operator

  Operator(Op_UC;
					 nr_orb=size(Op_UC,1), 
					 acts_on_atoms=atom,
					 kwargs...)

end

# ----- inverse participation ratio ----- #


function IPR(;kwargs...)::Function

	ldos = LDOS(;kwargs...)

	return function ipr(P::AbstractMatrix{<:Number}; kw...)::Matrix{Float64}

		sum(abs2, ldos(P; kw...), dims=ldos.csdim) .|> inv 

	end 

	# special function, cannot be constructed directly from Operator(...)
	# = 1 if psi localized on a single site, 
	# = nr_atoms for psi spread on all atoms

end 


#=

function Velocity(BlochH, dir=1; dk=pi/100)

isnothing(dir) | isnothing(BlochH) && return (P;k) -> zeros(size(P,2),1)

function v(k) 

k1, k2 = copy(k), copy(k)

k1[dir] -= dk/2

k2[dir] += dk/2

return (BlochH(k2) - BlochH(k1))/dk

end

#	(BlochH(k .+ dk*(axes(k,1).==dir)) - BlochH(k))/dk

return (P;k,kwargs...) -> reshape(LA.diag(P'*v(k)*P),:,1)

end






function Partial_PositExpect_fromLDOS(
					Prob::AbstractMatrix, atoms::AbstractMatrix, str::AbstractString;
					kwargs...)

	Partial_PositExpect_fromLDOS(Prob, atoms,
															 Decode_PositionExpVal(str, err=true)...;
															 kwargs...)

end







function Partial_PositExpect_fromLDOS(
								Prob::AbstractMatrix, atoms::AbstractMatrix,
								D::Int, f::Function=identity; 
								convolute=false, delta=0.3, vsdim=2,
								kwargs...)

# P[energy, atom] = LDOS at energy, on atom 
# or: P[wf_index, atom] = localiz.prob on atom of a certain wf

#	A = the coordinates 'D'. B =  the perpendicular



	D_perp = ifelse(D==1, 2, ifelse(D==2, 1, nothing)) 


	function out(get_ylim_Pf, param_B)

		ylim, Prob_and_fA = get_ylim_Pf()

		exp_vals = zeros(length(param_B), size(Prob,1)) 
	
		for i in axes(exp_vals,1)

			Pi,fi = Prob_and_fA(i)

			@warn "WFs on columns" 

			exp_vals[i,:] = ArrayOps.Normalize(Pi,1; dim=dim)*fi

		end

		#	out[energy,some_y] = 	(sum over i in atoms_at_y)(P[energy, i] * f(x[i]))
		#												/(sum ...) P[energy, i]
	

#		exp_vals = ArrayOps.Normalize_Columns(exp_vals, 1)



		return (param_B, exp_vals, ["x","y","z"][D_perp], extrema(ylim))

	end


	sel(X,i) = selectdim(X, [2,1][dim], i)

	if !convolute


		uniq_B, inds_B = Utils.Unique(sel(atoms, D_perp),
																	inds="all", sorted=true)

		return out(uniq_B) do 

			f_A = f.(sel(atoms, D))

			return (f_A, function(i) 

						inds_B[i] |> I -> (sel(Prob, I), f_A[I])

									end)

		end


	end  



	


	(A,dense_A,w_A), (B,dense_B,w_B) = map((D, D_perp)) do d
	
		v = sel(atoms, d)
	
		u = Utils.Unique(v, sorted=true)

		m,M = extrema(u)

#		m,M = [m,M] .+ (d==dim)*[1,-1]*1/4*(M-m) # disregard edges
	
		return (v, range(m, M, length=150), minimum(diff(u))*delta)
	
	end


	return out(dense_B) do 

		convProb = Prob*Algebra.get_CombinedDistribution(
												(A, B),
												(repeat(dense_A, outer=length(dense_B)),
																 repeat(dense_B, inner=length(dense_A))),
												(w_A, w_B);
												weights="Lorentzian", normalize=false)

		LI = LinearIndices((axes(dense_A,dim), axes(dense_B,dim)))

		f_A = f.(dense_A)

		return (f_A, function(i)
					s = LI[CartesianIndices((axes(dense_A,dim), i))][:]

					(sel(convProb,s), f_A)

				end)
	end

#	return out(get_ylim_Pf, dense_B)




#	L = Lorentzian with withd dB/5 or dA/5
#	P(R) = mapreduce(+, axes(atoms,1)) do i 
#				mapreduce(a->L(a[1]-a[2]), *, zip(R,atoms[i,:]), init=LDOS[i])
#			end

			
end








function Norm(;kwargs...)

	Operator([1];kwargs...)

end


function Overlap(;kwargs...)

  Operator([1],purpose="matrixelement";kwargs...)

end


#===========================================================================#
#
# Process input and construct the operator with a specific structure
# 
#---------------------------------------------------------------------------#


#===========================================================================#
#
# Operator for the case O = 1_{atoms}⊗ 1_{orbitals}
# 
#---------------------------------------------------------------------------#





#===========================================================================#
#
# Repeat an orbital-operator for several atoms
#	or an atom-operator for several orbitals
#
#---------------------------------------------------------------------------#



#===========================================================================#
#
# Normalize wavefunctions on columns
# 
#---------------------------------------------------------------------------#







=#






#############################################################################
end
