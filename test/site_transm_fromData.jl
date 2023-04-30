import myLibs: GreensFunctions, LayeredSystem, Lattices, Algebra, Utils, H_Superconductor,TBmodel, Operators, ObservablesFromGF, BandStructure, ArrayOps, Graph
import PyPlot, JLD ,SparseArrays , Combinatorics 
using LinearAlgebra, SharedArrays

include("mock_DevLeads.jl") 

colors = hcat(["brown","red","coral","peru"],["gold","olive","forestgreen","lightseagreen"],["dodgerblue","midnightblue","darkviolet","deeppink"]) 


colors = [colors[i,j] for i in axes(colors,1) for j in [3:size(colors,2);1:2]]

observables = [#"QP-Caroli", #"QP-BondVectorTransmission",
								"QP-SiteVectorTransmission",]
PyPlot.close.(1:10)


function get_hoppings(Bonds, Hopping::AbstractDict, 
											Slicer::Function, VirtLeads::AbstractDict=Dict(),
											)::Vector{Matrix{ComplexF64}}


	h = LayeredSystem.get_hoppings(Hopping[:Hopping], Slicer, VirtLeads)

	return Matrix{ComplexF64}[h(iR...)' for iR in zip(Bonds...)]

#											[Hopping[:Hopping](Rj,Ri) for (Ri,Rj) in Bonds[2]]

end 
 	
#hoppings = LayeredSystem.get_hoppings(dev_Hopping[:Hopping], 
#																			Slicer, 
#																			VirtLeads)
#
#
#BondHoppings = [hoppings(iR...)' for iR in zip(inds_DevBonds, Rs_DevBonds)]
println()

NR_ORB = rand(1:10)
	LENGTH = rand(1:30)
	WIDTH = rand(1:30)

#	NR_ORB = 1 
#	WIDTH=2
#	LENGTH=3 
#
	@show NR_ORB LENGTH WIDTH
	
	
	hopp = Dict(:Hopping=>get_hopp(NR_ORB), :nr_orb=>NR_ORB)
	
	
	DevLatt = get_latt([LENGTH, WIDTH])
	
	DevAtoms = Lattices.PosAtoms(DevLatt)


	inds_bonds_old = [C.I for C=findall(LinearAlgebra.triu!(isBond(DevAtoms,DevAtoms),1))]


	nr_bonds = length(inds_bonds_old) 


	Rs_bonds_old = Algebra.bondRs_fromInds(inds_bonds_old, DevAtoms; dim=2)

	Bonds_old = (inds_bonds_old,Rs_bonds_old) 




	
	leads = get_two_leads(max(1,div(WIDTH,2)), DevAtoms; hopp...)
	
	println()
	
	
for default_LAR in ["forced",
										default_LayerAtomRels(DevLatt)]
	
	NewGeom = LayeredSystem.NewGeometry(
									DevAtoms, default_LAR; 
									isBond=isBond, dim=2, nr_orb=NR_ORB,
									Leads=leads)
	
	LayerAtom,Slicer,LeadRels,VirtLeads=NewGeom 
	
	BondHoppings_old = get_hoppings(Bonds_old, hopp, Slicer, VirtLeads) 


	ial = LayerAtom[:IndsAtomsOfLayer]

	data_H,data_B= LayeredSystem.condStore_sharedHB(LayerAtom, DevAtoms; hopp...)

	data_A = LayeredSystem.prep_sharedA(LayerAtom, DevAtoms)


	PH = vcat(ones(div(NR_ORB,2)),zeros(NR_ORB-div(NR_ORB,2))) 
	
	projectors = [
		Operators.Projector(PH, :atoms; nr_orb=NR_ORB, dim=2),
		Operators.Projector(1 .- PH, :atoms; nr_orb=NR_ORB, dim=2),
		Operators.Projector(ones(NR_ORB), :atoms; nr_orb=NR_ORB, dim=2),
		]


println()
@testset "same bonds" begin 
	@test nr_bonds==LayeredSystem.nr_bonds(data_B)

	for ((l1,l2),B_new) in data_B, (i,j) in eachcol(B_new) 

		I = ial(l1)[i]
		J = ial(l2)[j]
		
		@test xor((I,J) in inds_bonds_old,(J,I) in inds_bonds_old)
	
	end 
end 
println() 


println()
@testset "same Hamilt" begin 

	for ((l1,l2),B_new) in data_B, (i,j) in eachcol(B_new) 

		I = ial(l1)[i]
		J = ial(l2)[j] 

		@test I!=J

		@test size(data_H[(l1,l2)])==(NR_ORB*length(ial(l1)),NR_ORB*length(ial(l2)))

		@test 1==count(zip(inds_bonds_old,BondHoppings_old)) do (b,h) 

			@test size(h)==(NR_ORB,NR_ORB)  

			@test LinearAlgebra.checksquare(h)==NR_ORB 

			if b==(I,J) #|| b==(J,I) || return false  

				return h≈TBmodel.slice_atoms(data_H[(l1,l2)],NR_ORB,i,j)'

			elseif b==(J,I)

				return h≈TBmodel.slice_atoms(data_H[(l1,l2)],NR_ORB,i,j)

			else 
				return false 
			end 

		end 

	end 
end 

	
	g_noH = LayeredSystem.LayeredSystem_toGraph(LayerAtom[:NrLayers], VirtLeads) 

	@show Graph.get_prop(g_noH,:UCsLeads)

	leadlabels = LayeredSystem.get_leadlabels(g_noH)
				
	
	
	E1 = rand() 
	
	
	
	G_old = GreensFunctions.GF_Decimation(hopp, VirtLeads, Slicer;
																							 LayerAtom...,
																							 leads_have_imag=true
																							 )(E1)
					
	G_new = GreensFunctions.GF_Decimation_fromGraph(E1, g_noH, data_H, Slicer;
																													leads_have_imag=true,
																													)

println()
@testset "G_new==G_old" begin 

	for i=1:LayerAtom[:NrLayers]

	@test G_new("Layer",1,"Layer",i)≈G_old("Layer",1,"Layer",i)

end 

end 
println()

	se = GreensFunctions.SelfEn_fromGDecim(G_new, VirtLeads, Slicer) 
	

	se(leadlabels[1],1)

	lead_ID = (leadlabels[1], 1)
									

	# data_H[(I,J)]: hopping from layer I to J  -- corresponding h(ri,rj)
	# SiteTransmission expects h(rj,ri)

	siteT_old = ObservablesFromGF.SiteTransmission(
										G_old,													# 1 arg 
										BondHoppings_old, Bonds_old..., # 3 args 
										se,															# 1 arg 
										lead_ID...;											# pos 6+7
										dim=1,
										)

	@show size(siteT_old)

	siteT_new = ObservablesFromGF.SiteTransmission(
										G_new,													# 1 arg 
										data_H, data_B, data_A..., # 4 args 
										se,															# 1 arg 
										lead_ID...;											# pos 7-8
										dim=1,
										nr_orb=NR_ORB,
										) 
println()
@testset "site_new == site_old" begin 
	@test size(siteT_new)==size(siteT_old) 

	@test LinearAlgebra.norm(siteT_new-siteT_old)<1e-8
#	@test siteT_new≈siteT_old 

end 
println() 



@testset "net T old==new" begin  

	for ps_=Combinatorics.powerset(projectors,0,2), ps=[ps_, reverse(ps_)]


tAB_new = ObservablesFromGF.CaroliConductance(G_new, se, leadlabels..., 1)
tAB_old = ObservablesFromGF.CaroliConductance(G_old, se, leadlabels..., 1)

@test tAB_new≈tAB_old   


end 
end 

for uc=1:Graph.get_prop(g_noH,:UCsLeads)-1

	ObservablesFromGF.CaroliConductance(G_new, se, leadlabels..., uc)

end 

end 

















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
	
	
	
	
	
	
	






































































































nothing
	
	
