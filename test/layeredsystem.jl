using myLibs: Lattices, LayeredSystem, Algebra,GreensFunctions, Graph,ArrayOps,TBmodel,Utils 

using LinearAlgebra, SparseArrays 
import PyPlot


function get_hopp(nr_orb_)

	intercell = rand(ComplexF64,nr_orb_,nr_orb_)
#(ones(nr_orb_,nr_orb_)-ArrayOps.UnitMatrix(nr_orb_))
	function hopp(ri,rj)

		local_pot = ArrayOps.UnitMatrix(nr_orb_)*isapprox(norm(ri-rj),0,atol=1e-8)*2.0  
	
		atom_hopp = intercell*isapprox(norm(ri-rj), 1, atol=1e-8) 


		return local_pot + atom_hopp 
	
	end  

end 

function get_latt(Nxy::Vector{Int})

	l1 = Lattices.SquareLattice()

	Lattices.Superlattice!(l1, Nxy, Labels=first) 

end 

function get_atoms(Nxy::Vector{Int})

	l2 = get_latt(Nxy)

	return Lattices.PosAtoms(l2)  

end 

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

function default_LayerAtomRels(DL::Lattices.Lattice)

	atoms = Lattices.PosAtoms(DL)
	layers = parse.(Int, Lattices.labelToComponent(DL)) .+ 1

		# layer of each atom [1, 2, ..., nr_at] => [L(1), L(2), ..., L(nr_at)]
	
	la,ial = Utils.FindPartners(enumerate(layers), sortfirst=true)

	return Dict(

			:NrLayers=> maximum(layers),

			:LayerOfAtom => la,

			:IndsAtomsOfLayer => ial,

			:AtomsOfLayer => function AtomsOfLayer(L::Int)::AbstractMatrix
			
									Lattices.Vecs(atoms, ial(L))

										end
							)
end




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


@testset "layered system hopping" begin 

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








end 




































nothing 
