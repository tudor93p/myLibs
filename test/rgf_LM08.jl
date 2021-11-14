include("rgf.jl")


PyPlot.close(8)



fig,Ax = PyPlot.subplots(2,2,num=8,figsize=(8,7))

set_title(fig)



for (i_term,term_x) in enumerate([:armchair,:zigzag])

	latt = ribbon_x(term_x, 3)

	for n=0:1
	
		for sl in Lattices.sublatt_labels(latt)


			plot_atoms(Ax[i_term,1], Lattices.Atoms_ManyUCs(latt; label=sl, Ns=n),label="$sl, n=$n")

		end 

	end 

end




[ax.legend() for ax in Ax[1:2,1]]


for (i_term,term_x) in enumerate([:armchair,:zigzag])

	LAR = slices_ribbon(term_x, nr_layers, nr_atoms_layer)

	plot_layers(Ax[i_term,2]; LAR...) 

	plot_atoms(Ax[i_term,2], hcat(LAR[:AtomsOfLayer](0),
																LAR[:AtomsOfLayer](nr_layers+1));
						 label="+/-1",c="orange")

end 
























