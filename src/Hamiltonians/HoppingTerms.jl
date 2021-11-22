module HoppingTerms
#############################################################################

import Base: <, ==
import ..LA


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


struct HamiltBasis 

	use_spin::Bool 

	use_Nambu::Bool

#	spin::Vector{Int}

#	charge::Vector{Int}  

	matrix_dim::Int


	function HamiltBasis(spin::Bool, Nambu::Bool, args...)::HamiltBasis
	
		HamiltBasis(Val(spin), Val(Nambu), args...)
	
	end 
	
	
	function HamiltBasis(spin::Val{false}, Nambu::Val{false}, args...
											)::HamiltBasis
	
		new(false, false, #[1], [1], 
				1)
	
	end 
	
	
	function HamiltBasis(spin::Val{false}, Nambu::Val{true}, 
#											 spins::AbstractVector{Int}=[1,-1],
											 args...
											)::HamiltBasis
	
#		@assert length(spins)==2
	
		new(false, true, #spins, [1,-1], 
				2)
	
	end 
	
	
	function HamiltBasis(spin::Val{true}, Nambu::Val{false}, 
											 args...)::HamiltBasis
	
		new(true, false, #[1,-1], [1,1], 
				2)
	
	end 
	
	
	function HamiltBasis(spin::Val{true}, Nambu::Val{true}, 
											 args...)::HamiltBasis
	
		new(true, true, #repeat([1,-1], outer=2), repeat([1,-1],inner=2), 
				4)
	
	end 

end 





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
	
	return HamiltBasis(any(use_spin, spin), any(use_Nambu, Nambu))

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

	if spin1==spin2 && Nambu2
		
		return Nambu1 ? identity : Nambu_doubling

	elseif !spin1 && spin2 && !Nambu1 
		
		return Nambu2 ? spinNambu_doubling : spin_doubling 

	end

	error("Enlarge target basis")

end














struct HoppingTerm 

	basis::HamiltBasis 

	get_hopp_f::Function # give parameters => f(ri,rj)

	nr_uc::Int # max reach of the hopping

	# info about what f(ri,rj) returns 
	
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

































































































#############################################################################
end 

