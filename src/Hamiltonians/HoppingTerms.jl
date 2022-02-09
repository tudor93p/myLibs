module HoppingTerms
#############################################################################

import Base: <, ==, zero
import ..LA, ..TBmodel, ..Utils

export HamiltBasis, HoppingTerm

#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


struct HamiltBasis 

	use_spin::Bool 

	use_Nambu::Bool

	spin::Vector{Int}

	charge::Vector{Int}  

	matrix_dim::Int


	function HamiltBasis(spin::Bool, Nambu::Bool, args...)::HamiltBasis
	
		HamiltBasis(Val(spin), Val(Nambu), args...)
	
	end 
	
	
	function HamiltBasis(spin::Val{false}, Nambu::Val{false}, args...
											)::HamiltBasis
	
		new(false, false, [1], [1], 1)
	
	end 
	
	
	function HamiltBasis(spin::Val{false}, Nambu::Val{true}, 
											 spins::AbstractVector{Int}=[1,-1],
											 args...
											)::HamiltBasis
	
#		@assert length(spins)==2
	
		new(false, true, spins, [1,-1], 2)
	
	end 
	
	
	function HamiltBasis(spin::Val{true}, Nambu::Val{false}, 
											 args...)::HamiltBasis
	
		new(true, false, [1,-1], [1,1], 2)
	
	end 
	
	
	function HamiltBasis(spin::Val{true}, Nambu::Val{true}, 
											 args...)::HamiltBasis
	
		new(true, true, repeat([1,-1], outer=2), repeat([1,-1],inner=2), 4)

# equivalent to kron(Nambu, spin)
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

function ==(b1::HamiltBasis, b2::HamiltBasis)::Bool 

	for k in propertynames(b1)
		
		getproperty(b1, k)==getproperty(b2, k) || return false

	end 

	return true 

end
		


function <(b::HamiltBasis, B::HamiltBasis)::Bool 

	for k in (:use_spin ,:use_Nambu)

		getproperty(b, k) && !getproperty(B, k) && return false 

	end 

	return b.matrix_dim < B.matrix_dim 

end 

function zero(basis::HamiltBasis)::Matrix{ComplexF64}

	zeros(ComplexF64, basis.matrix_dim, basis.matrix_dim)

end 


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function zeroHoppTerm(basis::HamiltBasis)::Function 

	TBmodel.constHoppTerm(zero(basis))
	
end 


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


# a hopping term is characterized by: use_spin, use_Nambu 
# the desired output spin and Nambu must be larger

function compatible_HamiltBasis(basis::HamiltBasis, args...)::HamiltBasis

	compatible_HamiltBasis(basis.use_spin, basis.use_Nambu, args...)

end




function compatible_HamiltBasis(use_sN::NTuple{2,Bool}, args...)::HamiltBasis

	compatible_HamiltBasis(use_sN..., args...)

end



function compatible_HamiltBasis(use_spin::Bool, use_Nambu::Bool,
																spin::Bool, Nambu::Bool
																)::HamiltBasis

	use_spin && @assert spin

	use_Nambu && @assert Nambu
	
	return HamiltBasis(use_spin|spin, use_Nambu|Nambu)

end



function compatible_HamiltBasis(use_spin::Bool, use_Nambu::Bool,
																spin::Symbol, args...)::HamiltBasis 

	if spin==:min 
		
		return compatible_HamiltBasis(use_spin, use_Nambu, use_spin, args...)

	elseif spin==:max 

		return compatible_HamiltBasis(use_spin, use_Nambu, true, args...)

	end 

end 


function compatible_HamiltBasis(use_spin::Bool, use_Nambu::Bool,
																spin::Bool, Nambu::Symbol)::HamiltBasis

	if Nambu==:min 

		return compatible_HamiltBasis(use_spin, use_Nambu, spin, use_Nambu)

	elseif Nambu==:max

		return compatible_HamiltBasis(use_spin, use_Nambu, spin, true)

	end 

end





#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#



function spin_doubling(h::Union{Number,AbstractMatrix{<:Number}}
											 )::Matrix{ComplexF64}

	out = zeros(ComplexF64, 2, 2) 

	out[1,1] = out[2,2] = only(h)

	return out 

end 

function Nambu_doubling(h::Number)::Matrix{ComplexF64}

	out = zeros(ComplexF64, 2, 2)

	out[1,1] = h 

	out[2,2] = -conj(h)

	return out 

end 

function Nambu_doubling(h::AbstractMatrix{<:Number})::Matrix{ComplexF64}

	n = LA.checksquare(h)

	out = zeros(ComplexF64, 2n, 2n)

	out[1:n,1:n] = h 

	out[n+1:end,n+1:end] = -conj(h)

	return out  

end 


function spinNambu_doubling(h::Union{Number,AbstractMatrix{<:Number}}
											 )::Matrix{ComplexF64}

	h1 = only(h)

	out = zeros(ComplexF64, 4, 4)

	out[1,1] = out[2,2] = only(h)

	out[3,3] = out[4,4] = -conj(only(h))

	return out 

end 




#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function upgrade(hb1::HamiltBasis, hb2::HamiltBasis)::Function

	upgrade(hb1.use_spin, hb1.use_Nambu, hb2.use_spin, hb2.use_Nambu)

end


function upgrade(spin1::Bool, Nambu1::Bool,
								 spin2::Bool, Nambu2::Bool)::Function

	if spin1==spin2 
	
		Nambu1==Nambu2 && return identity 

		Nambu2 && !Nambu1 && return Nambu_doubling 

	elseif !spin1 && spin2 && !Nambu1 
		
		return Nambu2 ? spinNambu_doubling : spin_doubling 

	end

	@show spin1 spin2 Nambu1 Nambu2 

	error("Enlarge target basis")

end








#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function get_basis_signature(basis::HamiltBasis; kwargs...)::Vector{Int}

	#	out = map([(1,1), (1,-1), (-1,1), (-1,-1)]) do (C,S) 
	
	map(zip(basis.charge, basis.spin)) do (C,S)

		mapreduce(*, (
									(:charge,		C),
									(:electron,	C>0),
									(:hole,			C<0),

									(:spin,			S),
									(:spinup,		S>0),
									(:spindown, S<0),

									)
							) do (k,v)

			use = get(kwargs, k, false) 
		
			isa(use,Bool) && use && return v 
			
			isa(use,Int) && use==1 && return v 

			return 1 

		end 

	end 

end 



function basis_info(basis::HamiltBasis)::Function 

	kw0 = Dict(
						"Charge" => :charge,

						"QP" =>  nothing,

						"E" => :electron,

						"H" => :hole, 

#						"Sz" => :spin

						)

	function get_vector(label::Union{AbstractChar,Symbol};
											kwargs...)::Vector{Int}

		get_vector(string(label); kwargs...)

	end 

	function get_vector(label::AbstractString; kwargs...)::Vector{Int}

		@assert haskey(kw0,label) "Label '$label' not understood" 

		K = kw0[label]
		
		kw = isnothing(K) ? () : Dict(K=>true)

		out = get_basis_signature(basis; kw...)
		
		return length(unique(out))==1 ? out[1:1] : out 

	end 

	return get_vector 

end 















#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#










struct HoppingTerm 

	basis::HamiltBasis 

	tij::Function # tij(ri, rj, parameters...)

	nr_uc::Int # max reach of the hopping

	cond::Union{Bool,Function}
	# info about what f(ri,rj) returns 

end 


function HoppingTerm(basis::HamiltBasis, tij::Function, nr_uc::Int)::HoppingTerm

	HoppingTerm(basis, tij, nr_uc, true)

end

#function spinProjector(ht::HoppingTerm)::Matrix{Int}
#
#	LA.diagm(ht.spin_proj)
#
#end 

#function ==(ht1::HoppingTerm, ht2::HoppingTerm)::Bool
#
#	for k in propertynames(ht1)
#
#		getproperty(ht1, k)==getproperty(ht2, k) || return false 
#
#	end 
#u
#	return true 
#
#end


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function upgraded_tij(basis0::HamiltBasis, basis::HamiltBasis,
										 tij::Function)::Function

	upgrade(basis0, basis) ∘ tij
	
end 



function upgraded_tij(ht::HoppingTerm, basis::HamiltBasis,
										 tij::Function=ht.tij)::Function

	upgraded_tij(ht.basis, basis, tij)
	
end 


function upgraded_tij(ht::HoppingTerm, basis::HamiltBasis, params...)

	upgraded_tij(ht, basis, function tij(ri::AbstractVector{<:Real}, 
																			 rj::AbstractVector{<:Real}
																			 )::Union{Number,
																							 AbstractMatrix{<:Number}}

								 							ht.tij(ri, rj, params...)

													end)
end 


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#




#upgrade(implemented_spinNambu, desired_spinNambu)






#function upgrade(ht1::HoppingTerm, ht2::HoppingTerm)::Function 
#
#	ht1==ht2 && return identity 
#
#	@assert ht1<ht2 "Expand basis"
#
#
#	!ht1.use_spin && ht2.use_spin && error()
#
#	@assert !ht1.use_Nambu && ht2.use_Nambu 
#
##	spin_id = ArrayOps.UnitMatrix(2) 
#
#	z1 = zeros(ComplexF64, ht1.matrix_dim, ht1.matrix_dim)
#
#	return function upgrade_(m1::AbstractMatrix{<:Number})::Matrix{ComplexF64}
#
#		[m1 z1; z1 -conj(m1)]
#
#	end 
#
#
#end 
#
#



function init_hopp_term(HT::HoppingTerm,
												hopp_cutoff::Float64,
												 dist::Real,
												 dist_tol::Float64,
												 basis::HamiltBasis,
												 params...;
												 )::Tuple{Union{Bool,Function}, Function, Int, Real}

	if Utils.everything_is_null(params; atol=hopp_cutoff)

		return (false, zeroHoppTerm(basis), 0, 0)

	else 

		return (Utils.fSame(dist, dist_tol) ∘ -,

						upgraded_tij(HT, basis, params...), 

						HT.nr_uc, 

						dist) 

	end

end

#, params...)





























































































#############################################################################
end 

