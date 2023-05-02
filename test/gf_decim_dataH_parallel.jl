using Revise,Test 
import myLibs: GreensFunctions, LayeredSystem, Lattices, Operators, ObservablesFromGF
using LinearAlgebra , BenchmarkTools

import Profile 
include("mock_DevLeads.jl") 


println()


	
@testset "" begin 

	for LENGTH in [10,140,250]

		@show LENGTH 


	WIDTH =max(1,div(LENGTH,2))
	NR_ORB= 2

	hopp = Dict(:Hopping=>get_hopp(NR_ORB), :nr_orb=>NR_ORB)
	

	DevLatt = get_latt([LENGTH, WIDTH])
	
	DevAtoms = Lattices.PosAtoms(DevLatt)

	
	leads = get_two_leads(max(1,div(WIDTH,2)), DevAtoms; delta=false, hopp...)
	
	println()
	

	
#for default_LAR in ["forced", default_LayerAtomRels(DevLatt)]

default_LAR = default_LayerAtomRels(DevLatt)

	NewGeom = LayeredSystem.NewGeometry(
									DevAtoms, default_LAR; 
									isBond=isBond, dim=2, nr_orb=NR_ORB,
									Leads=leads)
	
	LayerAtom,Slicer,LeadRels,VirtLeads=NewGeom 
	

	data_single = LayeredSystem.condStore_sharedHB(LayerAtom, DevAtoms; hopp...)
	@time "single" data_single = LayeredSystem.condStore_sharedHB(LayerAtom, DevAtoms; hopp...)


#	data_multi = LayeredSystem.condStore_sharedHB(LayerAtom, DevAtoms; parallel=true, hopp...)
#	@time "multi" data_multi = LayeredSystem.condStore_sharedHB(LayerAtom, DevAtoms; parallel=true, hopp...)


#
#	PH = vcat(ones(div(NR_ORB,2)),zeros(NR_ORB-div(NR_ORB,2))) 
#	projectors = [ Operators.Projector(q, :atoms; nr_orb=NR_ORB, dim=2, sum_orbitals=false, sum_atoms=false) for q in [PH, 1 .- PH, 1 ] ]
#
#	data_A = LayeredSystem.prep_sharedA(LayerAtom, DevAtoms)
#
#
##	@show Base.summarysize(data_H)
##	@show Base.summarysize(data_B)
#
#	
#	g_noH = LayeredSystem.LayeredSystem_toGraph(LayerAtom[:NrLayers], VirtLeads) 
#
#
#	leadlabels = LayeredSystem.get_leadlabels(g_noH)
#				
#	
#	
#	E1 = rand() -0.5 + 0.002im 
#	
#
##	@show E1 
#
#	
#	G_ret = GreensFunctions.GF_Decimation_fromGraph(E1, g_noH, data_H, Slicer; leads_have_imag=false)
#	G_adv = GreensFunctions.GF_Decimation_fromGraph(conj(E1), g_noH, data_H, Slicer; leads_have_imag=false)
#
#println()
#
#println()
#
#	se = GreensFunctions.SelfEn_fromGDecim(G_ret, VirtLeads, Slicer) 
#
#
#	lead_ID = (leadlabels[1], 1)
##BenchmarkTools.DEFAULT_PARAMETERS.samples = 3 
##BenchmarkTools.DEFAULT_PARAMETERS.seconds = 60 
#	args = (G_ret,data_H, data_B, data_A..., se, lead_ID...) 
#	kwargs = (dim=1, nr_orb=NR_ORB)
#
##	@benchmark ObservablesFromGF.SiteTransmission($(args)...; $(kwargs)...)
#	ObservablesFromGF.SiteTransmission(args...; kwargs...)
#	


	
#	
#	for proj in projectors 
#
#		@time "DOS"	ObservablesFromGF.ComputeDOSLDOS_Decimation(G_new, LayerAtom[:NrLayers], LayerAtom[:IndsAtomsOfLayer], proj,  VirtLeads; dos=true, ldos=true, dim=2)
#	
#	
#		end  
#
#
end  end 



#@time f(5)
#@time f(6)
#@time f(7) 
#
#Profile.clear_malloc_data()
#f(140)
nothing 

