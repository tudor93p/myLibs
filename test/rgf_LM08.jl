include("rgf.jl")


PyPlot.close(8)

function ribbon_x_armchair(NrAtomsPerLayer::Int)

	make_y_ribbon!(honeycomb_x_armchair(), NrAtomsPerLayer)

end 





fig,axes = PyPlot.subplots(1,3,num=8,figsize=(12,4))

set_title(fig)


L0 = make_y_ribbon!(Lattices.HoneycombLattice(),2)

for n=1:4

	plot_atoms(axes[1], Lattices.Atoms_ManyUCs(L0;Ns=n), label="n=$n")

	plot_atoms(axes[2], Lattices.Atoms_ManyUCs(ribbon_x_armchair(3),Ns=n),label="n=$n")

end

axes[1].legend() 
axes[2].legend() 


for layer in 1:nr_layers

	atoms = LayerAtomRels[:AtomsOfLayer](layer)

	color = in(layer,1:nr_layers) ? (isodd(layer) ? "black" : "gray") : "red"

	plot_atoms(axes[3], atoms; color=color)

end 






























#
#println()#"\r                     ")
#
#PyPlot.figure(4)
#
#
#
