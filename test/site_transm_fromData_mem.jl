import myLibs: GreensFunctions, LayeredSystem, Lattices, Operators, ObservablesFromGF
using LinearAlgebra 

#include("mock_DevLeads.jl") 


println()

@testset "mem" begin 

	NR_ORB = rand(1:10)
	LENGTH = rand(2:30)
	WIDTH = rand(2:30)

	NR_ORB = 2 
	LENGTH = 100##200
	WIDTH = div(LENGTH,2)

	@show NR_ORB LENGTH WIDTH
	
	
	hopp = Dict(:Hopping=>get_hopp(NR_ORB), :nr_orb=>NR_ORB)
	
	
	DevLatt = get_latt([LENGTH, WIDTH])
	
	DevAtoms = Lattices.PosAtoms(DevLatt)

	
	leads = get_two_leads(max(1,div(WIDTH,2)), DevAtoms; hopp...)
	
	println()
	
	PH = vcat(ones(div(NR_ORB,2)),zeros(NR_ORB-div(NR_ORB,2))) 

	projectors = [ Operators.Projector(q, :atoms; nr_orb=NR_ORB, dim=2, sum_orbitals=false, sum_atoms=false) for q in [PH, 1 .- PH, 1 ] 
		]
	
for default_LAR in ["forced",
										default_LayerAtomRels(DevLatt)]
	
	NewGeom = LayeredSystem.NewGeometry(
									DevAtoms, default_LAR; 
									isBond=isBond, dim=2, nr_orb=NR_ORB,
									Leads=leads)
	
	LayerAtom,Slicer,LeadRels,VirtLeads=NewGeom 
	

	data_H,data_B= LayeredSystem.condStore_sharedHB(LayerAtom, DevAtoms; hopp...)

	data_A = LayeredSystem.prep_sharedA(LayerAtom, DevAtoms)



	
	g_noH = LayeredSystem.LayeredSystem_toGraph(LayerAtom[:NrLayers], VirtLeads) 


	leadlabels = LayeredSystem.get_leadlabels(g_noH)
				
	
	
	E1 = rand() 
	
	
	
	G_new = GreensFunctions.GF_Decimation_fromGraph(E1, g_noH, data_H, Slicer;
																													leads_have_imag=true,
																													)

println()

	for i=1:LayerAtom[:NrLayers]

	@test G_new("Layer",1,"Layer",i) isa AbstractMatrix



end 
println()

	se = GreensFunctions.SelfEn_fromGDecim(G_new, VirtLeads, Slicer) 
	

	se(leadlabels[1],1)

	lead_ID = (leadlabels[1], 1)
									

	siteT_new = ObservablesFromGF.SiteTransmission(
										G_new,													# 1 arg 
										data_H, data_B, data_A..., # 4 args 
										se,															# 1 arg 
										lead_ID...;											# pos 7-8
										dim=1,
										nr_orb=NR_ORB,
										)  

	@show size(siteT_new)
	
	println() 
	
	
	for proj in projectors 
	
		for dos=(true,false),ldos=(true,false)
	
			data = ObservablesFromGF.ComputeDOSLDOS_Decimation(G_new, LayerAtom[:NrLayers], LayerAtom[:IndsAtomsOfLayer], proj,  VirtLeads; dos=dos, ldos=ldos, dim=2)
	
	#		dos && @show data[1]
	#		ldos && @show extrema(data[2])
	
		end 
	
	
		end 
	end 
	nothing 
end 	
