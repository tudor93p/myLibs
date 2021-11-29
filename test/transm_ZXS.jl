# PHYSICAL REVIEW B 95, 245433 (2017)
# Quantum perfect crossed Andreev reflection in top-gated quantum anomalous
# Hall insulator–superconductor junctions
# Ying-Tao Zhang,1 Zhe Hou,2 X. C. Xie,2,3 and Qing-Feng Sun2,3,*

import myLibs: Algebra, H_Superconductor,  HoppingTerms, Utils, TBmodel
import myLibs.HoppingTerms: HoppingTerm, HamiltBasis
import LinearAlgebra


ZXS_inter = HoppingTerm(HamiltBasis(true,false),
						function inter_hopping(ri::AbstractVector,rj::AbstractVector,
																	 aux
																	 )::Matrix
						
							R = ri-rj 
						
							d = only(findall(isapprox.(abs.(R),1,atol=1e-5)))
						
							t = -Algebra.PauliMatrix(3) - 0.5im*Algebra.PauliMatrix(d)
						
							return R[d] > 0 ? t : t'
						
						 end, 1) 

ZXS_intra = HoppingTerm(HamiltBasis(true,false),

						function get_intra_hopping(ri::AbstractVector, rj::AbstractVector, m::Real, mu_L::Real)
						
							(m + 4)*Algebra.PauliMatrix(3)+ mu_L*Algebra.PauliMatrix(0)
						
						 end, 0)



function ZXS_hamilt(Delta, m, mu)

	hopp_cutoff = 1e-6 
	dist_tol = 1e-5 

	target_basis, sc_gap = H_Superconductor.init_SC_gapfunction(H_Superconductor.swave, Delta; spin_basis=true) 


  Append!, Sum, Hoppings = TBmodel.Add_Hopping_Terms(target_basis.matrix_dim,
																										 hopp_cutoff)
		# value, condition, (function )
  
  maxdist, d_update! = TBmodel.Update_Maxdist!(dist_tol, hopp_cutoff) 

  nmax, n_update! = TBmodel.Update_Maxdist!(0, hopp_cutoff)






	cond, val, n, d = HoppingTerms.init_hopp_term(ZXS_inter, hopp_cutoff, 1, dist_tol, target_basis, 1)

		Append!(Hoppings, 1, cond, val)

		n_update!(nmax, cond isa Function, n) 

		d_update!(maxdist, cond isa Function, d) 


	lp = HoppingTerms.init_hopp_term(ZXS_intra, hopp_cutoff, 
											0, dist_tol, target_basis, 
											m, mu
											)

	Append!(Hoppings, 1, lp[1], lp[2])





		val, fun, n = sc_gap

		Append!(Hoppings, val, true, fun)

		n_update!(nmax, val, n)

		d_update!(maxdist, val, n)
  
		

		return Dict(
		:Hopping => Sum(Hoppings,<(only(maxdist))∘LinearAlgebra.norm∘-),
		:Nr_UCs => Int(round(only(nmax))),	# provided the kind of hoppings above 
		:hopp_cutoff => hopp_cutoff,
		:dist_tol => dist_tol,
		:SC_basis => target_basis.use_Nambu,
		:BasisInfo => HoppingTerms.basis_info(target_basis),
		:nr_orb => target_basis.matrix_dim,
		)

end



f = ZXS_hamilt(0,0,0)[:Hopping]

for r in ([0,0],[1,0],[0,1],[2,1])


	@show r
	
	println.(eachrow(f([0,0],r)))

	println()

end




#2a: mu=1.5; Delta=?;

# m: min_m,max_m = [-1,1]*sqrt(Delta^2 + mu^2)












