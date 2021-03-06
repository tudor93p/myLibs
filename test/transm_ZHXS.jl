# PHYSICAL REVIEW B 95, 245433 (2017)
# Quantum perfect crossed Andreev reflection in top-gated quantum anomalous
# Hall insulatorâ€“superconductor junctions
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
		:Hopping => Sum(Hoppings,<(only(maxdist))âˆ˜LinearAlgebra.normâˆ˜-),
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

prep_lead = GreensFunctions.PrepareLead

#prep_lead(label, latt, lead_hopp::AbstractDict, del::Vararg{Float64}) = prep_lead(label, latt, lead_hopp, lead_hopp, del...)
#
#function prep_lead(label, latt::Lattices.Lattice, lead_hopp::AbstractDict, coupl_hopp::AbstractDict, del::Vararg{Float64})
#
#	LayeredSystem.PrepareLead(label, 
#														latt, 
#														Utils.add_args_kwargs(TBmodel.HoppingMatrix;
#																						 coupl_hopp...),
#														Utils.add_args_kwargs(TBmodel.HoppingMatrix; 
#																						 lead_hopp...),
#														GreensFunctions.GF_Surface(latt, "+", lead_hopp, del...))
#
#end 
#
#



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




function GF(dev_hopp, (LayerAtom,Slicer,LeadRels,VirtLeads); kwargs...)

	GreensFunctions.GF_Decimation(dev_hopp, VirtLeads, Slicer;
																LayerAtom..., kwargs...)

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


function plot_arrows(ax, XY_, dXY_; min_abs_len=nothing,
										 min_rel_len=nothing,
										 inds=Colon(),kwargs...)

	XY = XY_[:,inds]

	dXY = dXY_[:,inds]



	norms = LinearAlgebra.norm.(eachcol(dXY))

	max_norm = if !isnothing(min_abs_len) 
		
									min_abs_len 

						elseif !isnothing(min_rel_len)
		
							maximum(norms)*min_rel_len

						else 

							-1

						end


	mask  = norms .> max_norm

	for (xy,dxy) in zip(eachcol(XY[:,mask]),eachcol(dXY[:,mask]))

		plot_arrow(ax, xy, dxy; kwargs...)

	end 

	sleep(0.005)

end 



function get_local_bonds(show_transversal::Bool,
												 show_longitudinal::Bool,
												 inds_bonds,
												 Rs_bonds; 
												 IndsAtomsOfLayer::Function,
												 kwargs...)::Function 
	
	function glb(layer::Int)::BitVector 

		i_slice = IndsAtomsOfLayer(layer) 

		i_prev_slice = if show_longitudinal 
		
									IndsAtomsOfLayer(layer-1)

								else 

									nothing 

								end
	
		return map(collect(1:length(inds_bonds))) do I
	
			Ra,Rb = Rs_bonds[I]
			i,j = inds_bonds[I] 
	
			if show_transversal && !show_longitudinal
				
				Ra[1]â‰ˆRb[1]	|| return false 
	
			elseif show_longitudinal && !show_transversal
	
				Ra[2]â‰ˆRb[2] || return false 
	
			end 
	
	
			for j0 in i_slice 
			
				j0==j || continue 
	
				if show_longitudinal  
					
					for i0 in i_prev_slice 
	
						i0==i && return true 
	
					end 
	
				end 
				
				if show_transversal
	
					for i0 in i_slice 

#						i0==i && @show (i,j)

						i0==i && return true 
	
					end 
	
				end 
	
			end 
	
			return false 
	
		end 
	end 

end



























function transversal_current(Gr::Function, Ga::Function,
														 gamma::Function,
														 layer::Int,
														 lead...;#lead::AbstractString,uc_lead::Int;
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
	
	ij = zeros(Int,2,n)


	for i_nu in 1:n

		a = inds[i_nu] 

		b = inds[i_nu+1] 

		ij[1,i_nu] = a 
		ij[2,i_nu] = b 

		Ra,Rb = selectdim(atoms, dim, i_nu), selectdim(atoms, dim, i_nu+1) 

		t_ab = Hopping(Ra, Rb)

		# t from nu==a to nu+1==b

		@assert any(abs.(t_ab).>0.1) t_ab

		ga_La = Ga(lead...,"Atom",a)

		gr_aL = Gr("Atom",a,lead...)

		gr_bL = Gr("Atom",b,lead...)

		ggg_ba = gr_bL*gamma(lead...)*ga_La


		ga_Lb = Ga(lead..., "Atom", b)


		@assert isapprox(ga_La,gr_aL',rtol=1e-8) && isapprox(ga_Lb,gr_bL',rtol=1e-8)



		ggg_ab = gr_aL*gamma(lead...)*ga_Lb

		cond = -1im*LinearAlgebra.tr(t_ab*ggg_ba - t_ab'*ggg_ab)  

		@assert isapprox(imag(cond),0,atol=1e-8) cond 


		selectdim(xy, dim, i_nu) .= Ra/2+Rb/2 

		selectdim(dxy, dim, i_nu) .= real(cond) * LinearAlgebra.normalize(Rb-Ra)

	end 



	return ij,xy,dxy 


end 



function longitudinal_current(Gr::Function, Ga::Function,
														 layer_::Int,
														 hopp::Function,
														 LeadLayerSlicer::Function,
														 proj1::Function=identity,
														 proj2::Function=identity;
														 AtomsOfLayer::Function,
														 IndsAtomsOfLayer::Function,
														 dim::Int,
														 kwargs...)

	dir_left,left_ = "left",-1
	current = 0 
	dir_right,right_ = "right",+1 

# equivalent to current due to the *left* lead 
	

#	dir_left,left_ = "right",1
#	current = 0
#	dir_right,right_ = "left",-1



	left,layer,right = [left_,current,right_] .+ layer_




	atoms_left = AtomsOfLayer(left)  

	atoms = AtomsOfLayer(layer)

	atoms_right = AtomsOfLayer(right) 


	sigma_left = GreensFunctions.SelfEn(hopp(atoms_left,atoms),
																			Gr("Layer",left,dir=dir_left)
																			)|>proj1
	
	sigma_right = GreensFunctions.SelfEn(hopp(atoms_right,atoms),
																			 Gr("Layer",right,dir=dir_right)
																			 )|>proj2


	gamma_right = GreensFunctions.DecayWidth(sigma_right)
#	gamma_left = GreensFunctions.DecayWidth(sigma_left)


	gr = Gr("Layer",layer) 

	ga = Ga("Layer", layer)

	
	
	ggg = gr*gamma_right*ga 

	@assert gr â‰ˆ ga' 

	
	cond_matrix = -im*(ggg*sigma_left' - sigma_left*ggg)



	inds = IndsAtomsOfLayer(layer)

	inds_left = IndsAtomsOfLayer(left)

	n = length(inds)
	
	@assert length(inds_left)==n

	xy = zeros(2,n)
	
	dxy = zeros(2,n)
	
	ij = zeros(Int,2,n)



	for (i_rel,(i_abs,i_abs_left)) in enumerate(zip(inds,inds_left))

		Ra = selectdim(atoms_left, dim, i_rel)
		Rb = selectdim(atoms, dim, i_rel)

		ij[:,i_rel] .= [i_abs_left,i_abs]


		L,slice = LeadLayerSlicer("Atom",i_abs) 

		L_l,slice_l = LeadLayerSlicer("Atom",i_abs_left)

		@assert L == ("Layer",layer)
		@assert L_l == ("Layer",left) 


		cond = LinearAlgebra.tr(cond_matrix[slice...,slice...])
	
		@assert isapprox(imag(cond),0,atol=1e-7) cond 



		selectdim(xy, dim, i_rel) .= Ra/2+Rb/2 

		selectdim(dxy, dim, i_rel) .= real(cond) * LinearAlgebra.normalize(Rb-Ra)

	end 

	return ij,xy, dxy

end 



best_linear_comb = Utils.best_linear_combs_int





