import myLibs: GreensFunctions, LayeredSystem, Lattices, Algebra, Utils, H_Superconductor,TBmodel, Operators, ObservablesFromGF

P1 = (length = 40, Barrier_height = 1.75, Barrier_width = 0.03, SCDW_phasediff = 0.78, SCDW_p = 2, SCDW_width = 0.005, SCDW_position = 0, SCpx_magnitude = 0.4, delta = 0.002, width = 20, SCpy_magnitude = 0.4, Hopping=1.0, ChemicalPotential=2.0)

P2 = (Attached_Leads = "AB",)

P3 = [(ChemPot = 0.0, Hopping = 1.0, Label = "A", Direction = -1, Contact = (-1, 1), Lead_Width = 10, Lead_Coupling = 1.0, SCbasis = true), (ChemPot = 0.0, Hopping = 1.0, Label = "B", Direction = 1, Contact = (1, -1), Lead_Width = 10, Lead_Coupling = 1.0, SCbasis = true)]




Dev = Lattices.SquareLattice()

NxNy = [P1[:length], P1[:width]]

Lattices.ShiftAtoms!(Dev, n=(1 .- NxNy)/2)

Lattices.Superlattice!(Dev, NxNy) 

Lattices.ReduceDim!(Dev)


DevAtoms = Lattices.PosAtoms(Dev) 

@show Lattices.NrVecs(DevAtoms)

SurfaceAtoms = Lattices.SurfaceAtoms(DevAtoms,1)

SurfaceAtoms2 = Lattices.filterAtoms_byNrBonds(-5:-1, DevAtoms, 1)


@assert isapprox(SurfaceAtoms,SurfaceAtoms2) 


@show Lattices.NrVecs(SurfaceAtoms)  

CornerAtoms = Lattices.filterAtoms_byNrBonds(-2, DevAtoms, 1)


@show Lattices.NrVecs(CornerAtoms)  

qs = Utils.Quadrant(CornerAtoms, Algebra.Mean(CornerAtoms, 2), dim=2) 

#@show qs 

CornerAtoms = CornerAtoms[:,sortperm(qs)]

#println.(eachcol(CornerAtoms))

LeadContacts = map(P3) do lead_params

#		2 ----------- 1
#		|      |      |     
#		|	---- : ----	|
#		|      |      |
#		3 ----------- 4


end 

LeadLatts = map(P3) do lead_params 

	(x,y) = lead_params[:Contact]


	t,(a,b) = if abs(x) == 1 
								
								(y,x>0 ? [4,1] : [3,2])

						elseif abs(y)==1 

								 (x,y>0 ? [2,1] : [3,4])

						end 

	start = (CornerAtoms[:,a]*(1-t) + CornerAtoms[:,b]*(1+t))/2

	@show start 


	d = lead_params[:Direction]::Int 

	sDir,Dir = sign(d),abs(d)

	L = Lattices.SquareLattice(string("Lead", sDir>0 ? "+" : "-", Dir))

	Lattices.Superlattice!(L, [1,lead_params[:Lead_Width]])


	Lattices.KeepDim!(L, Dir) 



	@show Lattices.NrVecs(Lattices.PosAtoms(L))
	@show Lattices.PosAtoms(L)[:,1]


	shift = start - Algebra.Mean(Lattices.PosAtoms(L), 2) 

	@show shift 
	
	Lattices.ShiftAtoms!(L, r=shift)

	Lattices.Align_toAtoms!(L, SurfaceAtoms)
	

	@show Lattices.PosAtoms(L)[:,1]

	return L

end 

println()

inds_DevBonds = Algebra.get_Bonds(DevAtoms, 1; dim=2) 

Rs_DevBonds =  Algebra.bondRs_fromInds(inds_DevBonds, DevAtoms; dim=2)

@show size(inds_DevBonds) size(Rs_DevBonds)

leadlabels = [string(l) for l in P2[:Attached_Leads]]

@show leadlabels


println()

function zerogap(args...; kwargs...)::Matrix{ComplexF64}

	zeros(ComplexF64, 4, 4)

end 


coupling_Hparam = map(P3) do lp
	
	(Hopping = lp[:Lead_Coupling], SC_Gap = (zerogap, 0),	)

end 

lead_Hparam = map(P3) do lp 
	
	(Hopping = lp[:Hopping], 
	 ChemicalPotential = lp[:ChemPot], 
	 SC_Gap = (zerogap, 0),
								)

end 




function chiral_gap(ri::AbstractVector, rj::AbstractVector)::Matrix{ComplexF64}

	R = ri/2 + rj/2 
	
	t = tanh( R[1] / (P1.SCDW_width * P1.length) )
	
	eta0 = [P1[:SCpx_magnitude], P1[:SCpy_magnitude]*1im]
	
	phase_factor = exp(0.5im*P1.SCDW_phasediff*pi*sign(t))
	
	eta_x,eta_y = eta0 .* setindex!(ones(2), t, P1.SCDW_p) .* phase_factor 


	delta(u)::Bool = isapprox(u,0,atol=0.0001)
	

	psi = 0  

	(x, y)= ri-rj 
	
	D = zerogap()
	
	D[1,4]=-1im*eta_y*(-delta(y - 1) + delta(y + 1))*delta(x)/2 + (1im*eta_x*delta(x - 1)/2 - 1im*eta_x*delta(x + 1)/2 + psi*delta(x))*delta(y)
	
	D[2,3]=-1im*eta_y*(-delta(y - 1) + delta(y + 1))*delta(x)/2 + (1im*eta_x*delta(x - 1)/2 - 1im*eta_x*delta(x + 1)/2 - psi*delta(x))*delta(y)
	
	D[3,2]=-1im*(-delta(y - 1) + delta(y + 1))*delta(x)*conj(eta_y)/2 + (1im*delta(x - 1)*conj(eta_x)/2 - 1im*delta(x + 1)*conj(eta_x)/2 - delta(x)*conj(psi))*delta(y)
	
	D[4,1]=-1im*(-delta(y - 1) + delta(y + 1))*delta(x)*conj(eta_y)/2 + (1im*delta(x - 1)*conj(eta_x)/2 - 1im*delta(x + 1)*conj(eta_x)/2 + delta(x)*conj(psi))*delta(y)
	
	return D
		
end 
			


dev_Hparam = (Hopping = P1[:Hopping],
							ChemicalPotential = P1[:ChemicalPotential],
							#LocalPotential = local_potential(dev_params),
							SC_Gap = (chiral_gap, 1))

dev_Hopping = H_Superconductor.SC_Domain(dev_Hparam,	[0,1])

@show size(dev_Hparam.SC_Gap[1]([10,0],[11,0])) 

@show size(lead_Hparam[1].SC_Gap[1]([10,0],[11,0])) 

leads = map(zip(leadlabels, LeadLatts, coupling_Hparam, lead_Hparam, P3)) do (
											label, latt, coupling_HP, lead_HP, lp)

	hopp_c, hopp_l = map([coupling_HP, lead_HP]) do Params 

		H_Superconductor.SC_Domain(Params,	[0,1])

	end 

	hc, hl = map([hopp_c, hopp_l]) do hopp

		function hm(args...) 
			
			TBmodel.HoppingMatrix(args...; hopp...)

		end 

	end 
	
	
	@assert hopp_l[:Nr_UCs]<=1 

  get_TBL(L::Lattices.Lattice) = Lattices.NearbyUCs(L, hopp_l[:Nr_UCs])

	intra, inter = TBmodel.BlochHamilt_ParallPerp(latt, get_TBL, hopp_l)

	# intra might be a function of k,
	# inter might be more than one matrix (m ∈ {1,2,...})

	function gf(E::Number; k=[])::Matrix{ComplexF64}

		GreensFunctions.GF_SanchoRubio(E, intra(k), only(values(inter));
																	 target="+")["+"]
		
	end

	return LayeredSystem.PrepareLead( label, latt, hc, hl, gf)

end 

observables = ["QP-BondVectorTransmission",
								"QP-SiteVectorTransmission",]

@show keys(leads[1])

isBond = Algebra.EuclDistEquals(1; tol=1e-5,dim=2)

NewGeom = LayeredSystem.NewGeometry(DevAtoms, "forced"; 
																		Leads=leads, isBond=isBond, dim=2,
																		dev_Hopping...)

@show length(NewGeom)


println()

LayerAtom,Slicer,LeadRels,VirtLeads=NewGeom 


G = GreensFunctions.GF_Decimation(
					dev_Hopping, VirtLeads, Slicer;
					LayerAtom...,
					)


hoppings = LayeredSystem.get_hoppings(dev_Hopping[:Hopping], Slicer, VirtLeads)


BondHoppings = [hoppings(iR...)' for iR in zip(inds_DevBonds, Rs_DevBonds)]

#G, Energy, observables, Slicer, VirtLeads 

trace_bond = Operators.Trace("atoms", [1,1,1,1]; dim=2, nr_orb=4,
																	sum_up=true, nr_at=1)



ENERGIES = LinRange(-1,1,100)|>collect 

out = Utils.RecursiveMerge(ENERGIES.+0.002im; dim=2) do Energy 

	print("\rEnergy=",real(round(Energy,digits=3)),"      ")

	g = G(Energy)

#	@show size(g("Atom",1,"Atom",2))

	self_energy = GreensFunctions.SelfEn_fromGDecim(g, VirtLeads, Slicer)

	return Dict(map(leadlabels) do lead

		bvt = ObservablesFromGF.BondTransmission(
																		i -> g("Atom", i, lead, 2),
																		self_energy(lead, 2),
																		inds_DevBonds, Rs_DevBonds, BondHoppings;
																		f=trace_bond) 

#out["QP-BondVectorTransmission"][lead] 

		return lead => ObservablesFromGF.SiteTransmission(bvt, inds_DevBonds, Rs_DevBonds)

	end)


end  

println("\r                             ")



@show keys(out) 

@show size.(values(out))

@show [maximum(abs, v) for v in values(out)]
