module HoppingTerms
#############################################################################

import Base: <, ==, zero
import ..LA, ..TBmodel, ..Utils, ..Algebra

export HamiltBasis, HoppingTerm

#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


const ORBITALS = [:spin, :Nambu]

const NR_ORBITALS = 2 

#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function orbital_matrix(spin::AbstractMatrix{<:Number},
												Nambu::AbstractMatrix{<:Number},
#												other_orbitals
												)::AbstractMatrix

	kron(Nambu, spin)

end 

function orbital_matrix(spin::AbstractVector{<:Number},
												Nambu::AbstractMatrix{<:Number}
												)::AbstractMatrix

	orbital_matrix(LA.Diagonal(spin), Nambu)

end 


function orbital_matrix(spin::AbstractMatrix{<:Number},
												Nambu::AbstractVector{<:Number}
												)::AbstractMatrix

	orbital_matrix(spin, LA.Diagonal(Nambu))
	
end  

function orbital_matrix(spin::AbstractVector{<:Number},
												Nambu::AbstractVector{<:Number}
												)::AbstractVector 

	LA.diag(orbital_matrix(LA.Diagonal(spin), LA.Diagonal(Nambu)))
	
end  





#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function order_orbitals(;kwargs...)::Tuple

	Tuple(kwargs[orb] for orb in ORBITALS)

end 

function order_orbitals(default_val; kwargs...)::Tuple

	Tuple(get(kwargs, orb, default_val) for orb in ORBITALS)

end 

#function order_orbitals(f::Function; kwargs...)
#
#	f((kwargs[orb] for orb in [:spin, :Nambu])...)
#
#end 
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




function orbital_signature(orbital::Symbol, used::Bool, active::Bool=true)::Vector

	used || return [1] 

	orbital==:spin && return active ? [1,-1] : [1,1]

	orbital==:Nambu && return active ? [1,-1] : [1,1]

	error("Orbital '$orbital' not implemented")

end  

function basis_size(orbital::Symbol, args...)::Int 

	length(orbital_signature(orbital, args...))

end 

function basis_size(status::Vararg{Bool,NR_ORBITALS})::Int
	
	prod(basis_size(ou...) for ou in zip(ORBITALS,status))

end   


function basis_size(;kwargs...)::Int 
	
	basis_size(order_orbitals(false; kwargs...)...)

end 



function orbital_signature(orbital::Symbol, used::Bool, active_orb::Symbol
													)::Vector

	orbital_signature(orbital, used, orbital==active_orb) 

end 



function orbital_signature_many(active_orb::Symbol,
																status::Vararg{Bool,NR_ORBITALS}
																)

	orbital_matrix((orbital_signature(ou..., active_orb) for ou in zip(ORBITALS,status))...)

end 

function orbital_signature_many(active_orbital::Symbol; orbitals...)

	orbital_signature_many(active_orbital,
												 order_orbitals(false; orbitals...)...)

end 


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


	function HamiltBasis(status::Vararg{Bool,NR_ORBITALS})::HamiltBasis

		new(status...,
				(orbital_signature_many(orb, status...) for orb in ORBITALS)...,
				basis_size(status...)
				)
	end 

	function HamiltBasis(;kwargs...)::HamiltBasis

		HamiltBasis(order_orbitals(false; kwargs...)...)

	end 


#	function HamiltBasis(spin::Bool, Nambu::Bool, args...)::HamiltBasis
#	
#		HamiltBasis(Val(spin), Val(Nambu), args...)
#	
#	end 
#	
#	
#	function HamiltBasis(spin::Val{false}, Nambu::Val{false}, args...
#											)::HamiltBasis
#
#
#		D = basis_size(;spin=false, Nambu=false)
#
#		s, n = [orbital_signature_many(orb; spin=false, Nambu=false) for orb in ORBITALS]
#			
#		@assert D==1 && s≈[1] && n≈[1] 
#
#		println("assert")
#
#		return new(false, false, [1], [1], 1)
#	
#	end 
#	
#	
#	function HamiltBasis(spin::Val{false}, Nambu::Val{true}, 
##											 spins::AbstractVector{Int}=[1,-1],
#											 args...
#											)::HamiltBasis
#	
##		@assert length(spins)==2
#
#		D = basis_size(;spin=false, Nambu=true)
#
#		spins, n = [orbital_signature_many(orb; spin=false, Nambu=true) for orb in ORBITALS]
#	
#		charges = [1,-1] # of the basis orbitals 
#
#		@assert D==2 && charges≈n
#
#
#		new(false, true, spins, charges , 2)
#	
#	end 
#	
#	
#	function HamiltBasis(spin::Val{true}, Nambu::Val{false}, 
#											 args...)::HamiltBasis
#	
#
#	spins = [1,-1]
#	charges = [1,1] 
#
#	D = basis_size(;spin=true)
#s, n = [orbital_signature_many(orb; spin=true) for orb in ORBITALS]
#	
#
#		@assert D==2 && charges≈n && spins≈s 
#
#		new(true, false, spins, charges, 2)
#	
#	end 
#	
#	
#	function HamiltBasis(spin::Val{true}, Nambu::Val{true}, 
#											 args...)::HamiltBasis
#
#		old_spin = repeat([1,-1], outer=2)
#		old_Nambu = repeat([1,-1],inner=2) 
#
#		spin_signature = orbital_matrix([1,-1], [1,1])
#		charge_signature = orbital_matrix([1,1],[1,-1])
#
#		D = basis_size(;spin=true,Nambu=true)
#
#s, n = [orbital_signature_many(orb; spin=true,Nambu=true) for orb in ORBITALS] 
#		@assert old_spin≈spin_signature≈s
#		@assert old_Nambu≈charge_signature ≈n
#		@assert length(spin_signature)==length(charge_signature)==4==D
#
#
#		new(true, true, spin_signature, charge_signature, 4)
#
#	end 

end 




#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function uses_orbital(basis::HamiltBasis, orbital::Symbol)::Bool 

	getproperty(basis, Symbol("use_$orbital"))

end 

function orbital_signature(orbital::Symbol, basis::HamiltBasis, args...
													 )::Vector 

	orbital_signature(orbital, uses_orbital(basis, orbital), args...) 

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



function representation(basis::HamiltBasis, args::Vararg{Any,NR_ORBITALS}
												)::AbstractVecOrMat

	matrices = map(zip(args, ORBITALS)) do (a,orb) 
		
		n = basis_size(orb, basis)
		
		a isa Number && return fill(a, n) 

		if a isa AbstractVecOrMat{<:Number}

			N = isa(a, AbstractVector) ? length(a) : LA.checksquare(a)

			N==1 && return fill(a[1], n)

			N==n && return a 

		end 

		return zeros(n)

	end 


	return orbital_matrix(matrices...)

end 

function representation(basis::HamiltBasis; kwargs...)::AbstractVecOrMat

	representation(basis, order_orbitals(hcat(1); kwargs...)...)

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

	function get_info(label::Union{AbstractChar,Symbol};
											kwargs...)::Vector{Int}

		get_info(string(label); kwargs...)

	end 

	function get_info(label::AbstractString; kwargs...)::Vector{Int}

		@assert haskey(kw0,label) "Label '$label' not understood" 

		K = kw0[label]
		
		kw = isnothing(K) ? () : Dict(K=>true)

		out = get_basis_signature(basis; kw...)
		
		return length(unique(out))==1 ? out[1:1] : out 

	end 

	function get_info(inds::Vararg{Any,NR_ORBITALS})::AbstractVecOrMat

		representation(basis, 
									 ((isa(i,Int) && 1<=i<=3) ? Algebra.PauliMatrix(i) : 1 for i in inds)...)

	end  



#
		
	function get_info(;kwargs...)::AbstractVecOrMat 
	
		get_info(order_orbitals(nothing; kwargs...)...)

	end 


	return get_info 

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


function HoppingTerm(basis::HamiltBasis, 
										 tij::Function, 
										nr_uc::Int=1)::HoppingTerm

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

function upgraded_tij(tij::Function)::Function 

	tij 
end 



function upgraded_tij(tij::Function,
											p1::Union{Number,Utils.List},
											params...
											)::Function 
	
	function tij_(ri::AbstractVector{<:Real},
								rj::AbstractVector{<:Real},
								)::Union{Number, AbstractMatrix{<:Number}}

				tij(ri, rj, p1, params...)

	end

end 



function upgraded_tij(basis0::HamiltBasis, basis::HamiltBasis,
										 tij::Function)::Function

	upgrade(basis0, basis) ∘ tij
	
end 



function upgraded_tij(ht::HoppingTerm, basis::HamiltBasis,
										 tij::Function=ht.tij)::Function

	upgraded_tij(ht.basis, basis, tij)
	
end 


function upgraded_tij(ht::HoppingTerm, basis::HamiltBasis, 
											p1::Union{Number,Utils.List},
											params...)::Function

	upgraded_tij(ht, basis, upgraded_tij(ht.tij, p1, params...))

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






























































































#############################################################################
end 

