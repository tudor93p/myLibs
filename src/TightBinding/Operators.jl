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

	data::Union{Number, <:VecOrMat{Float64}, <:VecOrMat{ComplexF64}, Function} 

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

		return new(data, diag, dim..., ws, get_iter_inds(ws, ao, numbers))

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

		@assert nr_at==nr_orb==1 "Confusing sizes. Provide 'diag' arg."

	end  

	return first(possib_diag)

end 


function prep_diag(A::Union{Number,<:AbstractVecOrMat{<:Number}},
									 numbers::NTuple{3,Int}, diag::Symbol)::Symbol

	possib_diag = which_diag(A, numbers) 

	@assert diag in possib_diag "The 'diag' given is not consistent with sizes"

	return diag 

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

	(dim, [2,1][dim])#, Utils.sel(dim), Utils.sel([2,1][dim]))

end 


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function get_acts_on(diag::Symbol, nr_at::Int, nr_orb::Int, size_H::Int;
										 kwargs...
#										 acts_on_atoms::AbstractVector{Int}=1:nr_at,
#										 acts_on_orbs::AbstractVector{Int}=1:nr_orb,
)::NTuple{2,Vector{Int}}

	(get(kwargs, :acts_on_atoms, 1:nr_at), 
	 get(kwargs, :acts_on_orbs, 1:nr_orb)
	 )

#  diag==:atoms && return 1:nr_at 
#
#	diag==:orbitals && return 1:nr_orb 
#
#
#  error("It's unclear on what the operator acts")

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
										::Val{:orbitals}
										)

	@assert only(unique(size(A)))==L2

end 

function check_size(A::AbstractVecOrMat{<:Number},
										::Tuple{Val{true}, Val{false}},
									 (L1,L2,L12)::NTuple{3,Int},
										::Val{:atoms}
										)

	@assert only(unique(size(A)))==L1

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

function check_size(A::Number,
										::NTuple{2,Val},
										arg,
										::Val{:all}
										)

end

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
									 args...)::Vector{<:Number}


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
									 args...)::Vector{<:Number}

	v = ExpectVal(P, Op, (Val(false),Val(false)), csdim, args...) 

	return dropdims(sum(v, dims=csdim), dims=csdim)

end  


function ExpectVal(P::AbstractMatrix{<:Number}, 
									 Op::AbstractMatrix{<:Number},
									 ::NTuple{2,Val{true}},
									 csdim::Int,
									 I::Colon,
									 args...)::Vector{<:Number}

	LA.diag(P' * Op * P)

end 




#===========================================================================#
#
# Build matrix when only one of atoms/orbitals is iterated 
#
#---------------------------------------------------------------------------#


function ExpectVal(P::AbstractMatrix{<:Number},
										iter_Op::Union{T, <:AbstractVector{T}},
									 csdim::Int,
									 I::AbstractVector{<:AbstractVector{Int}}
									 )::Matrix{<:Number} where T<:AbstractVecOrMat{<:Number}

	Utils.VecsToMat(iter_Op, I; dim=csdim) do Op,inds

		ExpectVal(P, Op, (Val(true),Val(true)), csdim, inds)

	end 

end   



function ExpectVal(P::AbstractMatrix{<:Number}, 
									 Op::Union{Number, <:AbstractVector{<:Number}},
									 sums::Union{Tuple{Val{true},Val{false}},
															 Tuple{Val{false}, Val{true}}},
									 csdim::Int,
									 I::AbstractVector{<:AbstractVector{Int}},
									 diag::Val,
									 args...)::Matrix{<:Number}

	ExpectVal(P, get_iter_Op(Op, sums, I, diag), csdim, I)

end  


function ExpectVal(P::AbstractMatrix{<:Number}, 
									 Op::AbstractMatrix{<:Number},
									 sums::Tuple{Val{true},Val{false}},
									 csdim::Int,
									 I::AbstractVector{<:AbstractVector{Int}},
									 diag::Val{:orbitals},
									 args...)::Matrix{<:Number}

	ExpectVal(P, get_iter_Op(Op, sums, I, diag), csdim, I)

end 



function ExpectVal(P::AbstractMatrix{<:Number}, 
									 Op::AbstractMatrix{<:Number},
									 sums::Tuple{Val{false},Val{true}},
									 csdim::Int,
									 I::AbstractVector{<:AbstractVector{Int}},
									 diag::Val{:atoms},
									 args...)::Matrix{<:Number}

	ExpectVal(P, get_iter_Op(Op,sums,I,diag), csdim, I)


end 






#===========================================================================#
#
# using an Operator as a one-argument function => compute expectation 
#
#---------------------------------------------------------------------------#

function (H::Operator)(P::AbstractMatrix)

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



function extend_data(A::T,
										 diag::Val,
										 numbers::NTuple{3,Int},
										 args...)::T where T<:AbstractVecOrMat{<:Number}

	check_size(A, (Val(false),Val(true)), numbers, diag) 

	return extend_data_(A, diag, numbers, args...)

end



function extend_data_(A::T, ::Val{:none}, args...
										 )::T where T<:AbstractVecOrMat{<:Number}
	A

end

function extend_data_(A::T,
										 ::Union{Val{:atoms},Val{:orbitals}},
										 ::NTuple{3,Int},
										 ::Union{Tuple{Val{true},Val{false}},
															 Tuple{Val{false},Val{true}}},
										 )::T where T<:AbstractVecOrMat{<:Number}
	A

end 


function extend_data_(A::T,
										 ::Val{:atoms},
										 numbers::NTuple{3,Int},
										 ::NTuple{2,V} where V<:Union{Val{true}, Val{false}}
										 )::T where T<:AbstractVecOrMat{<:Number}

	repeatOperator_manyAtoms(A, numbers...)

end




function extend_data_(A::T,
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

	# Op simulates a Hamiltonian with nr_at atoms and one orbital
	# 	will be reproduced for nr_orb orbitals

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
									kwargs...)

	a = selectdim(Rs, [2,1][dim], ax)

	Op = applicable(fpos, a) ? fpos(a) : fpos.(a)

	return Operator(Op, :orbitals; dim=dim, nr_at=nr_at, kwargs...)

end 





#===========================================================================#
#
# Trace over orbitals/atoms with some projection operator
# 
#---------------------------------------------------------------------------#



function Trace(what::AbstractString, Op=[1]; 
							 sum_up::Bool=false, kwargs...)::Function

	oa = ["orbitals","atoms"]

	what in oa || error("The first argument must be either of "*string(oa))


	S = size(Op,1)

	if sum_up && (S==1 || all(isapprox.(Op[1], Op[2:end], atol=1e-8))) 
		
		return A -> Op[1]*((typeof(A)<:AbstractMatrix) ? LA.tr : sum)(A)

	end


	(nr_at, nr_orb, size_H), enough_data = get_nratorbH(; kwargs...)


	!enough_data && return function (A)

								Trace(what, Op; sum_up=sum_up, size_H=size(A,1), kwargs...)(A)

									end

	Op = get_diag(Op)


	function check(A) 
		
		@assert size(A,1)==size_H 

		return A

	end


#	nr_at==1 && what=="atoms"





	S==size_H && sum_up && return A -> sum(Op.*get_diag(check(A)))


	what2 = oa[findfirst(!isequal(what), oa)]


	inds = TBmodel.Hamilt_indices_all(1:nr_orb,1:nr_at; iter=what2)


	iter(M) = get_diag(M) |> m -> (m[i] for i in inds)


	S==1 && !sum_up && return A -> Op[1]* sum.(iter(check(A)))


	S==size_H && return A -> [sum(oi.*ai) for (oi,ai) in zip(iter(Op),iter(check(A)))]





	sum_up_ = sum_up ? sum : identity


	nr = Dict("orbitals"=>nr_orb, "atoms"=>nr_at)



	S==nr[what] && return A -> sum_up_([sum(Op.*ai) for ai in iter(check(A))])

	S==nr[what2] && return A -> sum_up_(Op.*sum.(iter(check(A))))


	error("Something went wrong")



end




#=



	# ----- local charge operator ----- #


function Charge(Op_UC,atom;kwargs...)

  Operator(Op_UC,nr_orb=size(Op_UC,1),purpose="matrixelement",acts_on_atoms=vcat(atom);kwargs...)

end

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


	# ----- inverse participation ratio ----- #


function IPR(;kwargs...)

  ldos = LDOS(;kwargs...)

  return (P;kwargs...) -> 1 ./ sum(ldos(P).^2,dims=2)
		# special function, cannot be constructed from Operator(...)
		# = 1 if psi localized on a single site, 
		# = nr_atoms for psi spread on all atoms
end



	# ----- expectation of the position operator ----- #


function Decode_PositionExpVal(str;err=false)
	
	good(A...) = any(isequal(str), A)
									 
	for (d,A) in enumerate(["X","Y","Z"])
	
		good(A) && return d, identity
	
		good("|$A|") && return d, abs
	
		good(A*"2", "|$A|2") && return d, abs2 
		
		good(A*"3") && return d, x->x^3

		good("|$A|3") && return d, x->abs(x)^3

		good(A*"4", "|$A|4") && return d, x->x^4 
	
	end 

	err && error("'$str' not understood")

	return nothing,nothing

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

			exp_vals[i,:] = Algebra.Normalize(Pi,1; dim=dim)*fi

		end

		#	out[energy,some_y] = 	(sum over i in atoms_at_y)(P[energy, i] * f(x[i]))
		#												/(sum ...) P[energy, i]
	

#		exp_vals = Algebra.Normalize_Columns(exp_vals, 1)



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
