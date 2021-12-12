using myLibs:Lattices, LayeredSystem, Algebra,GreensFunctions


atoms = Lattices.PosAtoms(Lattices.Superlattice(Lattices.SquareLattice(), [10,7]))

isBond = Algebra.EuclDistEquals(1.0; dim=2) 


D, = LayeredSystem.LayerAtomRels(atoms, "forced"; isBond=isBond, dim=2) 


@show D[:NrLayers]
@show D[:LayerOfAtom](10)
@show D[:IndsAtomsOfLayer](3)
@show D[:AtomsOfLayer](3)




leads = Dict(:RightLead=>Dict(:label=>"A"), :LeftLead=>Dict(:label=>"B"))	 

g = LayeredSystem.LayeredSystem_toGraph(D[:NrLayers]; leads...)


#LayeredSystem.Plot_Graph(("test10",D[:IndsAtomsOfLayer]) ,D[:NrLayers],g)




println()
println()

Lead = Lattices.KeepDim!(Lattices.SquareLattice("Lead-A"), 1)

label = "A" 

@show Lattices.LattVec(Lead) |> size 

@show GreensFunctions.PrepareLead(label, Lead)[:head]


bridge = rand(2,1)

@show GreensFunctions.PrepareLead(label, Lead, bridge)[:head] .|>size


hopp(x,y=x) = rand(2,2)

@show GreensFunctions.PrepareLead(label, Lead, hopp, hopp, hopp)[:intracell]

@show GreensFunctions.PrepareLead(label, Lead, bridge, hopp, hopp, hopp)[:GF](rand()) .|>size


















