import myLibs: GreensFunctions, LayeredSystem, Lattices, Algebra, Utils, H_Superconductor,TBmodel, Operators, ObservablesFromGF, BandStructure, ArrayOps
import PyPlot, JLD ,SparseArrays , LinearAlgebra

include("mock_DevLeads.jl") 

colors = hcat(["brown","red","coral","peru"],["gold","olive","forestgreen","lightseagreen"],["dodgerblue","midnightblue","darkviolet","deeppink"]) 


colors = [colors[i,j] for i in axes(colors,1) for j in [3:size(colors,2);1:2]]

observables = [#"QP-Caroli", #"QP-BondVectorTransmission",
								"QP-SiteVectorTransmission",]
PyPlot.close.(1:10)




@testset "" begin 
	NR_ORB = 3 
	LENGTH = 20 
	WIDTH = 10 
	
	
	
	
	
	DevLatt = get_latt([LENGTH, WIDTH])
	
	DevAtoms = Lattices.PosAtoms(DevLatt)
	
	
	hopp = (Hopping=get_hopp(NR_ORB), nr_orb=NR_ORB)
	
	leads = get_two_leads(max(1,div(WIDTH,2)), DevAtoms; hopp...)
	
	println()
	
	
	
#	default_LAR = "forced"
	
	default_LAR = default_LayerAtomRels(DevLatt) 
	
	NewGeom = LayeredSystem.NewGeometry(
									DevAtoms, default_LAR; 
									isBond=isBond, dim=2, nr_orb=NR_ORB,
									Leads=leads)
	
	LayerAtom,Slicer,LeadRels,VirtLeads=NewGeom 
	
	data_H,data_B = LayeredSystem.condStore_sharedHB(LayerAtom, DevAtoms; hopp...)
	
	PH = vcat(ones(div(NR_ORB,2)),zeros(NR_ORB-div(NR_ORB,2)))
	el_proj = Operators.Projector(PH, :atoms; nr_orb=NR_ORB, dim=2)
	ho_proj = Operators.Projector(1 .- PH, :atoms; nr_orb=NR_ORB, dim=2)
	

	
	g_noH = LayeredSystem.LayeredSystem_toGraph(LayerAtom[:NrLayers], VirtLeads) 
	
	leadlabels = LayeredSystem.get_leadlabels(g_noH)
	
	
	E1 = rand() 
					
	G_new = GreensFunctions.GF_Decimation_fromGraph(E1, g_noH, data_H, Slicer;
																													leads_have_imag=true,
																													)

	G_new("Layer",1,"Layer",1)

	
	se = GreensFunctions.SelfEn_fromGDecim(G_new, VirtLeads, Slicer) 

	se(leadlabels[1],1)
	
	
	#function get_hoppings(Bonds, Hopping::AbstractDict, 
	#											Slicer::Function, VirtLeads::AbstractDict=Dict(),
	#											)::Vector{Matrix{ComplexF64}}
	#
	#
	#	h = LayeredSystem.get_hoppings(Hopping[:Hopping], Slicer, VirtLeads)
	#
	#	return [h(iR...)' for iR in zip(Bonds...)]
	#
	##											[Hopping[:Hopping](Rj,Ri) for (Ri,Rj) in Bonds[2]]
	#
	#end 
			
	#hoppings = LayeredSystem.get_hoppings(dev_Hopping[:Hopping], 
	#																			Slicer, 
	#																			VirtLeads)
	#
	#
	#BondHoppings = [hoppings(iR...)' for iR in zip(inds_DevBonds, Rs_DevBonds)]
	
		
	#	for i in 1:nr_at, j in i+1:nr_at
	#	
	#		I = TBmodel.Hamilt_indices(1:nr_orb, i, nr_orb)  
	#		J = TBmodel.Hamilt_indices(1:nr_orb, j, nr_orb)  
	#	
	#	
	#		Tij = ObservablesFromGF.BondTij(GWG[I,J], H_ll[J,I])
	#																		
	#		ObservablesFromGF.addSiteTiTj!(layer_siteT, 2, 
	#																	 (i,j), 
	#																	 (A0[:,i], A0[:,j]),
	#																	 Tij) 
	#	
	#		if i_atom_in_layer in [i,j]
	#	
	#			print_transm(inds_atoms_of_layer[[i,j]], (A0[:,i], A0[:,j]), Tij, out1)
	#	
	#		end 
	#	
	#	end  
	
	
		
		#gA(i)::AbstractMatrix{ComplexF64} = g_("Atom",i,"A",2)
	#	
	#	bvt = ObservablesFromGF.BondTransmission(g_, 
	#																					 BondHoppings,
	#																		inds_DevBonds, 
	#																		self_energy_, "A",2;
	##																		f=trace_bond) 
	#															)																		
	#	
	#	
	#	
	#	mask = [i_atom in I for I in inds_DevBonds]
	#	
	#	out2 = Dict() 
	#	
	#																		f=trace_bond, 
#																			dim=1)
		#@show size(jx0_)
		
	#	close_to_dw = isapprox.(DevAtoms[1,:], minimum(abs,DevAtoms[1,:]))  
		
		
	#	jx = ObservablesFromGF.SiteTransmission_fromBondT(bvt, inds_DevBonds, Rs_DevBonds; dim=1)
	#	
	#	println("Site transmission: ",round.(jx[i_atom,:],digits=4)) 
	#	
	#	println()
	#	PyPlot.close(4)
	#	error()  
	#	
	#	@show transpose(jx[LayerAtom[:IndsAtomsOfLayer](layer),:])
	#	
	#	
	#	jx_ = ObservablesFromGF.SiteTransmission(g_, 
	#																					 BondHoppings,
	#																		inds_DevBonds, Rs_DevBonds,
	#																		self_energy_, "A",2;
	#																		dim=1)
	##
	##	
	#
	#	error() 
	#	
	#	
	#	PyPlot.figure(4)
	#	
	#	#PyPlot.quiver(eachrow(DevAtoms)..., eachcol(jx)...; pivot="middle")
	#	
	#	PyPlot.quiver(eachrow(hcat(DevAtoms,A0))..., eachrow(hcat(transpose(jx),siteT))...; pivot="middle")
	#	
	#	
	#	
	#	
	#	gSD = g_("A", 2, "B", 2)
	#	
	#	
	#	sigmaS = GreensFunctions.SelfEn_fromGDecim(g_, VirtLeads, Slicer)("A", 2)
	#	sigmaD = GreensFunctions.SelfEn_fromGDecim(g_, VirtLeads, Slicer)("B", 2)
	#	
	#	cm = ObservablesFromGF.longitudinal_conductance_matrix(gSD, sigmaS, sigmaD)
	#	
	#	
	#	
	#	neighb_uc = uc#in [uc-1,uc+1]
	#	
	#	
	#	bondT = SparseArrays.spzeros(Float64, 
	#															 size(atoms(uc),2), 
	#															 size(atoms(neighb_uc),2))
	#	
	#	
	#	
	#	Gi = g_("A", uc, "A", neighb_uc)
	#	
	#	A = Gi*GreensFunctions.DecayWidth(self_energy_("A",neighb_uc))*Gi' 
	#	
	#	
	#	B = lA[:intracell][1]  
	#	
	#	#only i>j  should be considered 
	#	
	#	nr_orb = 4
	#	
	#	for j in axes(bondT,2)
	#	
	#		J = TBmodel.Hamilt_indices(1:nr_orb, j, nr_orb)  
	#	
	#		for i in axes(bondT,1)
	#	
	#			I = TBmodel.Hamilt_indices(1:nr_orb, i, nr_orb) 
	#	
	#			Tij = ObservablesFromGF.BondTij(view(A,I,J), view(B,J,I))
	#	
	#	#		transversal formula gives nonzero for lead intra-cell 
	#	#		transversal formula gives zeros for lead inter-cell 
	#	
	#	
	#	
	#			isapprox(Tij,0,atol=1e-14) || setindex!(bondT, Tij, i, j)
	#	
	#		end 
	#	
	#	end 
	#	
	#	@show SparseArrays.findnz(bondT)
	#	
	#	for i in 1:nr_at, j in i+1:nr_at 
	#	
	#		ObservablesFromGF.addSiteTiTj!(siteT, 2, (i,j), (A0[:,i],A0[:,j]), bondT[i,j])
	#	
	#	end 
	#	
	#	for (bond_index,(i,j)) in enumerate(lead_inter_bond_inds) 
	#	
	#		ObservablesFromGF.addSiteTiTj!(siteT, 2, (i,j), 
	#																	 (A0[:,i], AL[:,j]), bondTl[bond_index])
	#	
	#	
	#	end 
	#	
	#	
	#	
	#	#PyPlot.quiver(eachrow(A0)..., eachrow(siteT)...; pivot="middle")
	#	
	#	
		
	
	#jA = ObservablesFromGF.SiteTransmission(g_, 
	#																				 BondHoppings,
	#																	inds_DevBonds, Rs_DevBonds,
	#																	self_energy_, "A",2;
	#																	f=trace_bond, dim=1)
	#
	#
	##ObservablesFromGF.SiteTransmission_fromBondT(real(tr_orb(cm)),
	#
	#
	#PyPlot.figure(4) 
	#
	#PyPlot.quiver(eachrow(atoms1/2+atoms2/2)..., 
	#							eachrow(ArrayOps.multiply_elementwise(real(tr_orb(cm)),atoms2-atoms1))...;
	#							pivot="middle")
	#
	#
	#PyPlot.quiver(eachrow(atoms1/2+atoms0/2)..., 
	#							eachrow(ArrayOps.multiply_elementwise(cm2,atoms1-atoms0))...;
	#							pivot="middle")
	#
	##PyPlot.close.(1:10)
	
	
	
	
	
	#jx = jx[close_to_dw,1] 
	#
	#
	#jx0 = ObservablesFromGF.SiteTransmission_fromBondT(bvt0, inds_DevBonds, Rs_DevBonds; dim=1)
	#
	#
	#
	#@assert isapprox(jx0,jx0_) 
	#
	#jx0 = jx0[close_to_dw,1] 
	#
	#@show unique(DevAtoms[1,close_to_dw])
	#
	#@show size(jx) 
	#
	#
	#
	##@show ObservablesFromGF.TunnelingConductance_LeeFisher(g_, "A", 2)
	##@show ObservablesFromGF.TunnelingConductance_LeeFisher(g_, "B", 2)
	#
	#
	#
	#
	#
	#
	#
	#PyPlot.figure(5)
	#y = DevAtoms[2,close_to_dw]
	#
	#PyPlot.plot(y,jx, label="j", c=colors[1])
	#PyPlot.plot(y,jx0, label="j0", c=colors[2])
	#
	#PyPlot.xlabel("y")
	#PyPlot.ylabel("jx")
	#
	#PyPlot.legend()
	#	println()
	
	
	
	
	
	
	
end 






































































































nothing
	
	
