using myLibs: Lattices, LayeredSystem, Algebra,GreensFunctions, Graph
using LinearAlgebra 
import PyPlot


get_atoms(Nxy::Vector{Int}) = Lattices.PosAtoms(Lattices.Superlattice(Lattices.SquareLattice(), Nxy))  

get_leadlatt(n::Int=1, label...) = Lattices.KeepDim!(
																					Lattices.Superlattice( 
																					 Lattices.SquareLattice(label...),
																					 [1,n]
																					 )
																					 , 1)
get_lead(L::Lattices.Lattice) = GreensFunctions.PrepareLead(only(Lattices.sublatt_labels(L)), L)

get_lead(n::Int, label::String) = GreensFunctions.PrepareLead(label,get_leadlatt(n, label)) 

d_nn = norm(diff(get_atoms([2,1]),dims=2))


isBond = Algebra.EuclDistEquals(d_nn; dim=2)



@testset "LAR basisc" begin  
	atoms = get_atoms([10,7])



D, = LayeredSystem.LayerAtomRels(atoms, "forced"; isBond=isBond, dim=2) 

@show D[:NrLayers]
@show D[:LayerOfAtom](10)
@show D[:IndsAtomsOfLayer](3)
@show D[:AtomsOfLayer](3)





#LayeredSystem.LayeredSystem_toGraph(D[:NrLayers]; 
#																		RightLead=Dict(:label=>"A"), 
#																		LeftLead=Dict(:label=>"B"))
#

end 

#LayeredSystem.Plot_Graph(("test10",D[:IndsAtomsOfLayer]) ,D[:NrLayers],g)




println()
println()

@testset "lead basics" begin 
	Lead = get_leadlatt(1,"Lead-A")
	
	label = "A" 
	
	@test size(Lattices.LattVec(Lead))==(2,1)
	
	@test size.(GreensFunctions.PrepareLead(label, Lead)[:head])==[(2,1)]
	
	
	
	bridge = rand(2,3)
	
	@test size.(GreensFunctions.PrepareLead(label, Lead, bridge)[:head])==[(2,3),(2,1)]


	hopp(x,y=x) = rand(2,2)
	
	@test size.(GreensFunctions.PrepareLead(label, Lead, hopp, hopp, hopp)[:intracell])==[(2,2)]

	@test size.(GreensFunctions.PrepareLead(label, Lead, bridge, hopp, hopp, hopp)[:GF](rand())) == [(2,2),(2,2)]


println() 
end 



#@testset "" begin 

	Nxy = [10,10]

	lead_width = max(div(minimum(Nxy),2),1)

	atoms = get_atoms(Nxy)

leadA = get_lead(Lattices.Align_toAtoms!(get_leadlatt(lead_width,"A"),atoms,-1))
leadB = get_lead(Lattices.Align_toAtoms!(get_leadlatt(lead_width,"B"),atoms,1))

#	PyPlot.close()
#	PyPlot.scatter(Lattices.eachcomp(atoms)...,label="Device")
#	PyPlot.scatter(Lattices.eachcomp(only(leadA[:head]))...,label="A")
#	PyPlot.scatter(Lattices.eachcomp(only(leadB[:head]))...,label="B")
#
#	PyPlot.legend()




	LayerAtomRels,LeadRels,VirtLeads = LayeredSystem.NewGeometry(atoms, "forced";
																														 isBond=isBond,
																														 Leads=[leadA,leadB],
																														 dim=2)

	nr_layers = LayerAtomRels[:NrLayers]

@show nr_layers 

g = LayeredSystem.LayeredSystem_toGraph(nr_layers) 


LayeredSystem.add_leads!(g; VirtLeads...)










#end 




































nothing 
