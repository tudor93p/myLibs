include("rgf.jl")



function square_latt_lead(contact::Int, dir::Int, label::AbstractString;
													AtomsOfLayer::Function, kwargs...)

	atoms = AtomsOfLayer(contact) 


	a = dir * (if Lattices.NrVecs(atoms)>1 
							 
							 only(Lattices.Distances(atoms,1))

							else 1

							end)

	atoms[1,:] .+= a

	return Lattices.SquareLattice(a, label, atoms, nothing, 1), abs(a)

end  


function honeycomb_latt_lead(contact::Int, dir::Int, label::AbstractString;
														 AtomsOfLayer::Function, 
														 kwargs...)

	atoms = mapreduce(AtomsOfLayer, Lattices.Cat,
										contact+dir :dir: contact+4dir)

	@show size(atoms)

	L = ribbon_x_armchair2(Lattices.NrVecs(AtomsOfLayer(1)))

	a = only(Lattices.Distances(atoms,1))

	return Lattices.Lattice(dir*L.LattVec,atoms,
													nothing,
													L.LattDims), a
end 


function lead_contacts(contact::Int,
											 lead::Lattices.Lattice, a::Real; 
											 AtomsOfLayer::Function, IndsAtomsOfLayer::Function, 
											 kwargs...)


	isBond = Algebra.EuclDistEquals(a; tol=1e-5, dim=2)

	lc = LayeredSystem.get_LeadContacts(AtomsOfLayer(contact), latt, isBond)

	return IndsAtomsOfLayer(contact)[lc[1]]

end 



#										= Utils.zipmap(zip([-1,1],leadlabels,[1,NrLayers])
#																		 ) do (dir,label,contact)

#.|> collect



#honeycomb_LLC = Utils.zipmap(zip([-1,1],leadlabels,[1,NrLayers])
#																		 ) do (dir,label,contact)





fig,Ax = PyPlot.subplots(2,2,num=10,figsize=(12,12),sharex=true,sharey=true) 

nr_layers = 10
nr_atoms_layer = 10

LayerAtomRels = slices_armchair(nr_layers, nr_atoms_layer)




for layer in 1:nr_layers

	atoms = LayerAtomRels[:AtomsOfLayer](layer)

	kw = layer==1 ? Dict(:label=>"Device") : ()

	for ax in Ax 

		plot_atoms(ax, atoms; color=isodd(layer) ? "black" : "gray",
							 kw...)

	end

end 


for i=1:2,j=1:2 

	ax = Ax[i,j] 

	for args in [(1, -1, "Left", i),
									 	(LayerAtomRels[:NrLayers], 1, "Right", j)]

		lab, get_latt = [("square", square_latt_lead), 
															 ("honeycomb", honeycomb_latt_lead)
															 ][args[4]]

		latt,a = get_latt(args[1:3]...; LayerAtomRels...)

		for n=0:1

			kw = n==0 ? Dict(:label=>string(args[3]," $lab")) : Dict()

			plot_atoms(ax, Lattices.Atoms_ManyUCs(latt; Ns=n); kw...)
		
		end 

	end 

	ax.legend()

end 









