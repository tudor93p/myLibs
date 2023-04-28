using myLibs: Lattices, LayeredSystem, Algebra,GreensFunctions, Graph,ArrayOps,TBmodel
using LinearAlgebra, SparseArrays 
import PyPlot


function get_hopp(nr_orb_)

	function hopp(ri,rj)

		local_pot = ArrayOps.UnitMatrix(nr_orb_)*isapprox(norm(ri-rj),0,atol=1e-8)*2.0  
	
		atom_hopp = (ones(nr_orb_,nr_orb_)-ArrayOps.UnitMatrix(nr_orb_))*isapprox(norm(ri-rj), 1, atol=1e-8)*1.0
	
		return local_pot + atom_hopp 
	
	end  

end 

get_atoms(Nxy::Vector{Int}) = Lattices.PosAtoms(Lattices.Superlattice(Lattices.SquareLattice(), Nxy))  

get_leadlatt(n::Int=1, label...) = Lattices.KeepDim!(
																					Lattices.Superlattice( 
																					 Lattices.SquareLattice(label...),
																					 [1,n]
																					 )
																					 , 1)
get_lead(L::Lattices.Lattice, args...) = GreensFunctions.PrepareLead(only(Lattices.sublatt_labels(L)), L, args...)

get_lead(n::Int, label::String, args...) = GreensFunctions.PrepareLead(label,get_leadlatt(n, label), args...)

d_nn = norm(diff(get_atoms([2,1]),dims=2))


isBond = Algebra.EuclDistEquals(d_nn; dim=2)



#@testset "LAR basisc" begin  
	nr_orb = 2 
	atoms = get_atoms([10,7])



D, = LayeredSystem.LayerAtomRels(atoms, "forced"; isBond=isBond, dim=2) 

	nr_layers = D[:NrLayers]

	for i=axes(atoms,2)

		@test 1<=D[:LayerOfAtom](i)<=nr_layers

	end  

	for l=1:nr_layers 

		@test all(1 .<= D[:IndsAtomsOfLayer](l) .<= size(atoms,2))

		for j=1:nr_layers
			
			@test xor(l==j,isdisjoint(D[:IndsAtomsOfLayer](l),
																D[:IndsAtomsOfLayer](j)
																))
		end 

		D[:AtomsOfLayer](l)
	end 



	LayeredSystem.LayeredSystem_toGraph(nr_layers)


	hopp = Dict(:Hopping=>get_hopp(nr_orb), :nr_orb=>nr_orb) 

#	HoppMatr(args...) = TBmodel.HoppingMatrix(D[:AtomsOfLayer].(args)...; hopp...)

	g = LayeredSystem.LayeredSystem_toGraph(nr_layers, VirtLeads) # no H


	data_H,data_B = LayeredSystem.condStore_sharedHB(nr_layers,
																					 D[:IndsAtomsOfLayer],
																					 atoms;
																					 hopp...)


	@show keys(data_H)
	@show keys(data_B)  

#end 

#LayeredSystem.Plot_Graph(("test10",D[:IndsAtomsOfLayer]) ,D[:NrLayers],g)



println()
error()

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



@testset "layered system geometry" begin 

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



	g1 = LayeredSystem.LayeredSystem_toGraph(nr_layers, VirtLeads)

	@show nr_layers g1 

end 


#@testset "layered system hopping" begin 

	Nxy = [9,5]
	nr_orb = 2  




	atoms = get_atoms(Nxy)
	lead_width = max(div(minimum(Nxy),2),1)


	hopp = Dict(:Hopping=>get_hopp(nr_orb), :nr_orb=>nr_orb)

	leadA = get_lead(Lattices.Align_toAtoms!(get_leadlatt(lead_width,"A"),atoms,-1),hopp,0.001)

	leadB = get_lead(Lattices.Align_toAtoms!(get_leadlatt(lead_width,"B"),atoms,1),hopp,0.001)

	@show keys(leadA) keys(leadB)



	LayerAtomRels,Slicer,LeadRels,VirtLeads = LayeredSystem.NewGeometry(
					atoms, "forced"; isBond=isBond, Leads=[leadA,leadB], dim=2,
					nr_orb=nr_orb)


	nr_layers = LayerAtomRels[:NrLayers] 

	HoppMatr(args...) = TBmodel.HoppingMatrix(LayerAtomRels[:AtomsOfLayer].(args)...; hopp...)

	g2 = LayeredSystem.LayeredSystem_toGraph(HoppMatr, nr_layers, VirtLeads)

	for layer=1:nr_layers

		@show size(LayeredSystem.get_graphH(g2, "Layer", layer))
		
		layer==1 && continue 
		
		@show size(LayeredSystem.get_graphH(g2, "Layer", layer-1, "Layer", layer))


	end 








#end 




































nothing 
