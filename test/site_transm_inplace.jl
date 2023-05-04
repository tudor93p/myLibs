using Revise,Test
using myLibs: GreensFunctions, LayeredSystem, Lattices, Operators, ObservablesFromGF
using LinearAlgebra , BenchmarkTools 

#import Profile 
include("mock_DevLeads.jl") 

function get_projectors(nr_orb_)

	PH = vcat(ones(div(nr_orb_,2)),zeros(nr_orb_-div(nr_orb_,2))) 

	[ Operators.Projector(q, :atoms; nr_orb=nr_orb_, dim=2, sum_orbitals=false, sum_atoms=false) for q in [PH, 1 .- PH, 1 ]  
		]

end 


println()

NR_ORB= 2

hopp = Dict(:Hopping=>get_hopp(NR_ORB), :nr_orb=>NR_ORB)

	
@testset "" begin 
	
	LENGTH = 6 

	WIDTH =max(1,div(LENGTH,2))

	DevLatt = get_latt([LENGTH, WIDTH])
	
	DevAtoms = Lattices.PosAtoms(DevLatt)

	leads = get_two_leads(max(1,div(WIDTH,2)), DevAtoms; delta=false, hopp...)
	
	
	default_LAR = default_LayerAtomRels(DevLatt)

	LayerAtom,Slicer,LeadRels,VirtLeads = LayeredSystem.NewGeometry(
									DevAtoms, default_LAR; 
									isBond=isBond, dim=2, nr_orb=NR_ORB,
									Leads=leads)
	


	data_H,data_B = LayeredSystem.condStore_sharedHB(LayerAtom, DevAtoms; hopp...)
	data_A = LayeredSystem.prep_sharedA(LayerAtom, DevAtoms)

	
	g_noH = LayeredSystem.LayeredSystem_toGraph(LayerAtom[:NrLayers], VirtLeads) 

	leadlabels = LayeredSystem.get_leadlabels(g_noH)

	lead_ID = (leadlabels[1], 1) 


	E1 = rand() -0.5 + 0.002im 
	

#	@show E1 

	
	G_old = GreensFunctions.GF_Decimation(g_noH, E1, Slicer, data_H; leads_have_imag=false)


	self_en = GreensFunctions.SelfEn_fromGDecim(G_old, VirtLeads, Slicer) 

	pgf4 = GreensFunctions.prep_GF_fE(LayerAtom[:NrLayers], VirtLeads; leads_have_imag=false) 

	
	SE_new4 = GreensFunctions.SelfEn_Decimation(pgf4..., E1, Slicer, data_H)


	for k in leadlabels
		for uc=1:2

			@test self_en(k, uc) isa AbstractMatrix 
			@test self_en(k, uc) ≈ SE_new4(k, uc)

		end 
	end 

	@info "Self En OK"

	args = (G_old,data_H, data_B, data_A..., self_en, lead_ID...) 
	kwargs = (dim=1, nr_orb=NR_ORB)

#	@benchmark ObservablesFromGF.SiteTransmission($(args)...; $(kwargs)...)
	
	siteT1 = ObservablesFromGF.SiteTransmission(args...; kwargs...)

	siteT2 = zeros(ObservablesFromGF.size_SiteTransm(DevAtoms; dim=1))

	ObservablesFromGF.SiteTransmission!(siteT2, args...; kwargs...)

	@test siteT1 ≈siteT2 




	
	
	for proj in get_projectors(NR_ORB)

		data = ObservablesFromGF.ComputeDOSLDOS_Decimation(G_old, LayerAtom[:NrLayers], LayerAtom[:IndsAtomsOfLayer], proj,  VirtLeads; dos=true, ldos=true, dim=2)
		
		s = ObservablesFromGF.size_LDOS_Decimation(LayerAtom[:NrLayers], LayerAtom[:IndsAtomsOfLayer]; dim=2, VirtLeads...)

		ldos = zeros(s)

		ObservablesFromGF.LDOS_Decimation!(ldos, G_old, LayerAtom[:NrLayers], LayerAtom[:IndsAtomsOfLayer]; proj=proj, dim=2, VirtLeads...)

		@test data[2]≈ldos 

	
	
		end  


end 











nothing 
