# PHYSICAL REVIEW B 95, 245433 (2017)
# Quantum perfect crossed Andreev reflection in top-gated quantum anomalous
# Hall insulator–superconductor junctions
# Ying-Tao Zhang,1 Zhe Hou,2 X. C. Xie,2,3 and Qing-Feng Sun2,3,*

import myLibs: Algebra, H_Superconductor,  HoppingTerms, Utils, TBmodel,Lattices, BandStructure, Operators, LayeredSystem,GreensFunctions, ObservablesFromGF
import myLibs.HoppingTerms: HoppingTerm, HamiltBasis
import LinearAlgebra, PyPlot


ZHXS_inter = HoppingTerm(HamiltBasis(true,false),
						function inter_hopping(ri::AbstractVector,rj::AbstractVector,
																	 aux
																	 )::Matrix
						
							R = ri-rj 
						
							d = only(findall(isapprox.(abs.(R),1,atol=1e-5)))
						
							t = -Algebra.PauliMatrix(3) - 0.5im*Algebra.PauliMatrix(d)
						
							return R[d] > 0 ? t : t'
						
						 end, 1) 

ZHXS_intra = HoppingTerm(HamiltBasis(true,false),

						function get_intra_hopping(ri::AbstractVector, rj::AbstractVector, m::Real, mu_L::Real)
						
							(m + 4)*Algebra.PauliMatrix(3)+ mu_L*Algebra.PauliMatrix(0)
						
						 end, 0)



function ZHXS_hopping(Delta, m, mu)

	min_m,max_m = [-1,1]*sqrt(Delta^2 + mu^2)

	NrMZM = if m<min_m 

							2 

					elseif m<max_m 

						1

					else 

						0

					end 




	hopp_cutoff = 1e-6 
	dist_tol = 1e-5 

	target_basis, sc_gap = H_Superconductor.init_SC_gapfunction(H_Superconductor.swave, Delta; spin_basis=true, Nambu_basis=true)

  Append!, Sum, Hoppings = TBmodel.Add_Hopping_Terms(target_basis.matrix_dim,
																										 hopp_cutoff)
		# value, condition, (function )
  
  maxdist, d_update! = TBmodel.Update_Maxdist!(dist_tol, hopp_cutoff) 

  nmax, n_update! = TBmodel.Update_Maxdist!(0, hopp_cutoff)






	cond, val, n, d = HoppingTerms.init_hopp_term(ZHXS_inter, hopp_cutoff, 1, dist_tol, target_basis, 1)

		Append!(Hoppings, 0.5, cond, val)

		n_update!(nmax, cond isa Function, n) 

		d_update!(maxdist, cond isa Function, d) 


	lp = HoppingTerms.init_hopp_term(ZHXS_intra, hopp_cutoff, 
											0, dist_tol, target_basis, 
											m, mu
											)

	Append!(Hoppings, 0.5, lp[1], lp[2])





		val, fun, n = sc_gap

		Append!(Hoppings, val/2, true, fun)

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
		:NrMZM => NrMZM,
		)

end





function device_lattice(n::Int)

	latt = Lattices.SquareLattice()
	
	Lattices.Superlattice!(latt, [1,n]; recenter=true)

	Lattices.KeepDim!(latt, 1)

	return latt 

end 


function slices_ribbon(NrLayers::Int, NrAtomsPerLayer::Int)

	LayeredSystem.LayerAtomRels((device_lattice(NrAtomsPerLayer),
																NrLayers), "sublatt"; test=true)[1]

end

#Lattices.plot_atoms(latt/atoms, ax/n)
#
#LayeredSystem.plot_layers(ax/n)


function square_latt_lead(contact::Int, dir::Int, label::AbstractString;
													AtomsOfLayer::Function, 
													kwargs...)

	atoms = AtomsOfLayer(contact) 

	atoms[1,:] .+= dir 

	return Lattices.SquareLattice(dir, label, atoms, nothing, 1)

end   


two_lead_args(;NrLayers::Int, kwargs...) =  [(1, -1, "Left"), (NrLayers, 1, "Right")] 

prep_lead(label, latt, lead_hopp, del) = prep_lead(label, latt, lead_hopp, lead_hopp, del)

function prep_lead(label, latt, lead_hopp, coupl_hopp, del)

	LayeredSystem.PrepareLead(label, 
														latt, 
														Utils.add_args_kwargs(TBmodel.HoppingMatrix;
																						 coupl_hopp...),
														Utils.add_args_kwargs(TBmodel.HoppingMatrix; 
																						 lead_hopp...),
														GreensFunctions.GF_Surface(latt, "+", lead_hopp, del))

end 





function NewGeometry(LayerAtomRels, LeadContacts, leads, dev_hopp)

	LayeredSystem.NewGeometry(rand(0,0), 
															 LayerAtomRels; 
															 LeadContacts=LeadContacts,
															 Leads=leads, dim=2, 
															 dev_hopp...)


end


function lead_contacts(contact::Int,
											 lead::Lattices.Lattice; 
											 AtomsOfLayer::Function, IndsAtomsOfLayer::Function, 
											 kwargs...)

	isBond = Algebra.EuclDistEquals(1; tol=1e-5, dim=2)

	lc = LayeredSystem.get_LeadContacts(AtomsOfLayer(contact), lead, isBond)

	return IndsAtomsOfLayer(contact)[lc]

end 




function GF(dev_hopp, (LayerAtom,Slicer,LeadRels,VirtLeads))

	GreensFunctions.GF_Decimation(dev_hopp, VirtLeads, Slicer;
																			 LayerAtom...)

end 




function get_SE(g, (LayerAtom,Slicer,LeadRels,VirtLeads))

	GreensFunctions.SelfEn_fromGDecim(g, VirtLeads, Slicer)

end  




function get_Bonds(; kwargs...)

	inds_bonds = LayeredSystem.get_Bonds(1; kwargs...,  dim=2) 

	Rs_bonds = LayeredSystem.bondRs_fromInds(inds_bonds; kwargs..., dim=2)

	return inds_bonds, Rs_bonds

end  


function device_atoms(;AtomsOfLayer::Function, NrLayers::Int, 
											IndsAtomsOfLayer::Function, kwargs...)

	iaol = IndsAtomsOfLayer.(1:NrLayers)

	at = zeros(2, sum(length, iaol))

	for (layer,inds) in enumerate(iaol)
	
		for (i,a) in zip(inds,eachcol(AtomsOfLayer(layer)))

			at[:,i] = a 

		end 
	
	end

	return at 

end  



function plot_arrow(ax1, xy, dxy; kwargs...)

	ax1.arrow((xy-dxy/2)..., dxy...; length_includes_head=true, width=0.1, head_width=0.41, head_length=0.2, kwargs...)

end 



function transversal_current(Gr::Function, Ga::Function,
														 gamma::Function,
														 layer::Int,
														 lead::AbstractString;
														 Hopping::Function,
														 AtomsOfLayer::Function,
														 IndsAtomsOfLayer::Function,
														 dim::Int,
														 kwargs...)

	inds = IndsAtomsOfLayer(layer) 

	atoms = AtomsOfLayer(layer) 

	n = length(inds)-1 

	xy = zeros(2,n)
	
	dxy = zeros(2,n)
	


	for i_nu in 1:n

		a = inds[i_nu] 

		b = inds[i_nu+1]

		Ra,Rb = selectdim(atoms, dim, i_nu), selectdim(atoms, dim, i_nu+1) 

		t_ab = Hopping(Ra, Rb)

		# t from nu==a to nu+1==b

		@assert any(abs.(t_ab).>0.1) t_ab

		ga_La = Ga(lead,1,"Atom",a)

		gr_aL = Gr("Atom",a,lead,1) 

		gr_bL = Gr("Atom",b,lead,1)

		ggg_ba = gr_bL*gamma(lead,1)*ga_La


		ga_Lb = Ga(lead,1,"Atom",b)


		@assert isapprox(ga_La,gr_aL',rtol=1e-8) && isapprox(ga_Lb,gr_bL',rtol=1e-8)



		ggg_ab = gr_aL*gamma(lead,1)*ga_Lb

		cond = -1im*LinearAlgebra.tr(t_ab*ggg_ba - t_ab'*ggg_ab)  

		@assert isapprox(imag(cond),0,atol=1e-8) cond 


		selectdim(xy, dim, i_nu) .= Ra/2+Rb/2 

		selectdim(dxy, dim, i_nu) .= real(cond) * LinearAlgebra.normalize(Rb-Ra)

	end 



	return xy,dxy 


end 



#function longitudinal_current()
#
#	inds1 = IndsAtomsOfLayer(layer)
#
#	inds2 = IndsAtomsOfLayer(layer+1)
#
#
#	n = length(inds1)
#	
#	@assert 	length(inds2)==n
#
#	xy = zeros(2,n)
#	
#	dxy = zeros(2,n)
#
#
##	left = ("Layer",layer)
##	right = (
#
#	for (i,(i1,i2)) in enumerate(zip(inds1,inds2))
#
#		left = ("Atom",i1)
#
#		right = ("Atom", i2)
#
#		s_L = sigma(left...) 
#	
#		g_R = GreensFunctions.DecayWidth(sigma(right...))
#	
#		b = Gr(left..., right...)*g_R*Ga(right..., left...)
#	
#		
#		cond = -1im*LinearAlgebra.tr(b*s_L' - s_L*b)
#		
#		@assert isapprox(imag(cond),0,atol=1e-8) cond 
#
#
#		selectdim(xy, dim, i) .= Ra/2+Rb/2 
#
#		selectdim(dxy, dim, i) .= real(cond) * LinearAlgebra.normalize(Rb-Ra)
#
#	end 
#
#end 







