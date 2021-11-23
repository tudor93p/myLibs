module H_Superconductor

#############################################################################

import ..LA 

import ..Utils, ..TBmodel, ..Algebra

import ..HoppingTerms 

using ..HoppingTerms: HamiltBasis, HoppingTerm 

#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

include("chiral_pwave.jl")
include("swave.jl")
include("hopping_Peierls.jl")
include("local_potential.jl")

#step 1: ht_i = Number or init_hopp_term(some_specific_params)

#step 2: feed param_H = {ht_i} to SC_Domain 

#step 3: SC_Domain decides what the common basis is 

#step 4: SC_Domain constructs hoppings: Add!(ht(basis..., param...))



#SC_Gap = H_Superconductor.ChiralPwaveHopp(eta, basis_dim=:min)

#get_spin_basis(basis::HamiltBasis, spin::Symbol)::Bool 

#spin==:min && return 





function init_SC_gapfunction(SCHT::HoppingTerm,
												 param::Union{<:Function, <:Number,
																			 <:AbstractVector, Tuple},
												 hopp_cutoff::Float64,
												 dist_tol::Float64;
												 spin_basis=:min,
												 Nambu_basis=:min,
#												 latt_const::Float64=1.0,
												 )::Tuple{HamiltBasis, Tuple{Int, Function, Int}}



	if Utils.everything_is_null(param; atol=hopp_cutoff)

		basis = HoppingTerms.compatible_HamiltBasis(false, false, 
																								spin_basis, Nambu_basis)

		return basis, (0, HoppingTerms.zeroHoppTerm(basis), 0)

	else 

		basis = HoppingTerms.compatible_HamiltBasis(SCHT.basis, 
																								spin_basis, Nambu_basis)

		delta = Utils.fSame(0, dist_tol)

		function gap_hopping(ri::AbstractVector{<:Real}, 
												 rj::AbstractVector{<:Real}
											 )::Matrix{ComplexF64}

			Delta = zero(basis)

			SCHT.tij(ri, rj, param, Delta, basis, delta) # in-place Delta 

			return Delta 

		end 

		return basis, (1, gap_hopping, SCHT.nr_uc)

	end 	

end 



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#



function init_hopp_term(HT::HoppingTerm,
												hopp_cutoff::Float64,
												 dist::Real,
												 dist_tol::Float64,
												 basis::HamiltBasis,
												 params...;
												 )::Tuple{Union{Bool,Function}, Function, Int, Real}

	if Utils.everything_is_null(params; atol=hopp_cutoff)

		return (false, HoppingTerms.zeroHoppTerm(basis), 0, 0)

	else 

		return (Utils.fSame(dist, dist_tol) ∘ -,

						HoppingTerms.upgraded_tij(HT, basis, params...), 

						HT.nr_uc, 

						dist) 

	end

end




#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

#SC_gap_args(chiral_pwave_gapfunction, chiral_pwave_nmax, eta, hopp_cutoff,dist_tol) 





#chiral_pwave = (chiral_pwave_gapfunction!, 1, false)
#



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#






#===========================================================================#
#
# hopp
#
#---------------------------------------------------------------------------#

function ToBasis(Nambu::Bool=false, Z...)::Function 

	#  Nambu ? h->[h Zero; Zero -conj(h)] : h->h  
	
	ToBasis(Val(Nambu), Z...)


end

function ToBasis(::Val{false}, args...)::Function 

	TBmodel.matrixval 

end  

function ToBasis(::Val{true}, Zero::Tz=0.0im
								)::Function where Tz<:Union{Number, AbstractMatrix}

	function tobasis(h::Union{Number,AbstractMatrix})::Matrix{ComplexF64}

		[h Zero; Zero -conj(h)]
		
	end 
	
	tobasis(h::Function)::Function = tobasis ∘ h


	function tobasis(v::Union{Number, AbstractMatrix}, 
									 h::Union{Number,AbstractMatrix})::Matrix{ComplexF64}

		tobasis(TBmodel.matrixval(v,h))

	end 

	function tobasis(v::Union{Number, AbstractMatrix}, 
									 h::Function)::Function 

		tobasis∘TBmodel.matrixval(v, h)

	end 


	return tobasis 

	
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
#
#
#---------------------------------------------------------------------------#



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

										

#PeierlsHopp = HoppingTerm(PeierlsPhaseFactor, 1, false, false)
#
#ChiralPwaveHopp = HoppingTerm(chiral_pwave_gapfunction, 1, false, true)
#
#SwaveHopp = HoppingTerm(swave_gapfunction, 0, false, true)
#
#LocalPotHopp = HoppingTerm(args_local_potential, 0, false, false)

function get_basis_gap(::Nothing)::Tuple{Bool, HamiltBasis,Nothing}
	
	false, HamiltBasis(false,false), nothing 

end 


function get_basis_gap((basis,gap)::Tuple{HamiltBasis,Tuple}
											 )::Tuple{Bool, HamiltBasis, Tuple}
	
	true, basis, gap

end 







function SC_Domain(param_H_::NamedTuple, dist::AbstractVector{<:Real}; 
									 indomain=nothing, nr_uc::Int=1, 
									 dist_tol::Float64=1e-5, hopp_cutoff::Float64=1e-6,
									 )::Dict{Symbol,Any}


        # ---------------- initialize parameters --------------------- #      





  default_param_H = (
											Hopping			= 1.0,
											ChemicalPotential = 0.0,#*using_SC_gap,
											LocalPotential 		= 0.0,
										        SC_Gap			= nothing,
														Peierls = (0.0,:x),
										#	Electric_Field_x	= 0.0,
										#	Electric_Field_y	= 0.0,
										#	Electric_Field_z	= 0.0,
											Lattice_Imbalance	= 0.0,
										#	Rashba_SOC		= 0.0,
															)

  param_H = Utils.Combine_NamedTuples(param_H_, default_param_H)

	using_SC_gap, target_basis, SC_Gap = get_basis_gap(param_H[:SC_Gap])



#	@show default_param_H[:ChemicalPotential]
#	@show get(param_H_, :ChemicalPotential, "x")
#	@show param_H[:ChemicalPotential]

        # ------------------- needed for the hoppings ----------------- #      
  Append!, Sum, Hoppings = TBmodel.Add_Hopping_Terms(target_basis.matrix_dim,
																										 hopp_cutoff)
		# value, condition, (function )
  
  maxdist, d_update! = TBmodel.Update_Maxdist!(dist_tol, hopp_cutoff) 

  nmax, n_update! = TBmodel.Update_Maxdist!(0, hopp_cutoff)

        # -------------- nearest neighbor ------------------------ #      

	if length(dist)>=1 

		cond, val, n, d = init_hopp_term(
													hopping_Peierls, 
													hopp_cutoff, 
													dist[1], dist_tol, 
													target_basis,
													param_H[:Hopping],
													prep_PeierlsPhaseFactor(param_H[:Peierls]...)...
													) 

		Append!(Hoppings, 1, cond, val)

		n_update!(nmax, cond isa Function, n) 

		d_update!(maxdist, cond isa Function, d) 


  end  



        # -------------- atom bias (lattice imbalance) ----------- #      

#  Hoppings = Append(Hoppings, param_H[:Lattice_Imbalance],
#		(ri,rj) -> same(ri,rj),
#		(ri,rj) -> 2*ri[sl] - sum(values(slcoord)),
#							)

        # ----------------- local potential --------------- #      
				
	lp = init_hopp_term(local_potential, hopp_cutoff, 
											0, dist_tol, target_basis, 
											param_H[:LocalPotential], param_H[:ChemicalPotential]
											)

	Append!(Hoppings, 1, lp[1], lp[2])


        # -------------- superconducting pairing ----------------- #      

	if using_SC_gap

		val, fun, n = SC_Gap

		Append!(Hoppings, val, true, fun)

		n_update!(nmax, val, n)

		d_update!(maxdist, val, n)

  end
        # ------------------- electric field -------------------- #      

#  Hoppings = Append(Hoppings, param_H[:Electric_Field_x],
#		(ri,rj) -> same(ri,rj),
#		(ri,rj) -> ri[1],
#							)
#
#
#  Hoppings = Append(Hoppings, param_H[:Electric_Field_y],
#		(ri,rj) -> same(ri,rj),
#		(ri,rj) -> ri[2],
#							)
#
#
#  Hoppings = Append(Hoppings, param_H[:Electric_Field_z]*(dim==3),
#		(ri,rj) -> same(ri,rj),
#		(ri,rj) -> ri[3],
#							)
#
#        # -------------- anti-Haldane nnn imaginary hopping ------ #      
#
#
#  nn(ri) = findall(same.(Algebra.FlatOuterDist(hcat(ri[1:2]...),AllAtoms[:,1:2]),dist[1]))
#
#  chirality(ri,rj) = sign(Algebra.cross(eachcol(hcat(ri,rj) .- AllAtoms[intersect(nn(ri),nn(rj))[1],:] )...)[3])
#
#
#  Hoppings = Append(Hoppings, param_H[:AntiHaldane],
#		(ri,rj) -> same(ri[1:2].-rj[1:2],dist[2]),
#		(ri,rj) -> 1im/sqrt(27)*ri[sl]*chirality(ri,rj),
#							)
#
#  maxdist = update(maxdist,param_H[:AntiHaldane],dist[2])

        # -------------- Kane-Mele spin-orbit coupling ----------- #      

#    if same(la.norm(dif[:2]),dist_nnn)*same(abs(dif[2]),0.)*neighbors.size == 0
#  push!(hoppings, TBmodel.Hopping_Term( 
#		(ri,rj) -> isapprox(LA.norm(ri.-rj),dist_nnn,atol=dist_tol),
#		param_H[:KaneMele_SOC],
#		,
#		tol = hopp_cutoff))

        # -------------- Rashba spin-orbit coupling -------------- #      
#
#  push!(hoppings, TBmodel.Hopping_Term( 
#		(ri,rj) -> isapprox(LA.norm(ri.-rj),dist_nnn,atol=dist_tol),
#		param_H[:Rashba_SOC],
#		,
#		tol = hopp_cutoff))
#


        # -------------- compute max_dist and sum up ------------- #      

  function cond(ri::AbstractVector,rj::AbstractVector)::Bool 

		(isnothing(indomain) || indomain(ri,rj)) && LA.norm(ri-rj)<only(maxdist)

	end
  

	return Dict(
		:Hopping => Sum(Hoppings,cond),
		:Nr_UCs => Int(round(only(nmax))),	# provided the kind of hoppings above 
		:hopp_cutoff => hopp_cutoff,
		:dist_tol => dist_tol,
		:SC_basis => target_basis.use_Nambu,
		:BasisInfo => HoppingTerms.basis_info(target_basis),
		:nr_orb => target_basis.matrix_dim,
		)
			

end




#function ParticleHole_Flavours(Nambu)
#
#  samplehopp = rand(Complex{Float64},1+Nambu,1+Nambu) # add spin if Nambu
#
#  Zero,One = Base.zero(samplehopp),Base.one(samplehopp) 
#
#  return LA.diag(ToBasis(Nambu,Zero)(One)) # Nambu doubling
#
#end


#
#
#function ParticleHole_Operator(param_H,nr_at::Int64)
#
#  return TBmodel.Constant_Operator(nr_at,Charge_Operator_UC(param_H))
#end
#
#
#function Spin_Operator(param_H,atoms,axis=3)
#
#  if isnothing(param_H[:SC_Gap]) 
#
#    return TBmodel.Constant_Operator(size(atoms,1),hcat(0))
#
#  end
#
#  S = Algebra.PauliMatrices()[axis]
#
#  M = ToBasis(true,Base.zero(S))(S)
#
#  return TBmodel.Constant_Operator(size(atoms,1),M)
#
#end








###############################################################################

end

