include("rgf.jl")

nr_layers = 25

nr_atoms_layer = 50




fig,Ax = PyPlot.subplots(1,2,num=10,figsize=(12,8))#,sharex=true,sharey=true) 

LayerAtomRels = slices_armchair(nr_layers, nr_atoms_layer)


plot_layers(Ax[1]; LayerAtomRels...)

#
#
#for i=1:2,j=1:2 
#
#	ax = Ax[i,j] 
#
#	for args in [(1, -1, "Left", i),
#									 	(LayerAtomRels[:NrLayers], 1, "Right", j)]
#
#		lab, get_latt = [("square", square_latt_lead), 
#															 ("honeycomb", honeycomb_latt_lead)
#															 ][args[4]]
#
#		latt,a = get_latt(args[1:3]...; LayerAtomRels...)
#
#		for n=0:1
#
#			kw = n==0 ? Dict(:label=>string(args[3]," $lab")) : Dict()
#
#			plot_atoms(ax, Lattices.Atoms_ManyUCs(latt; Ns=n); kw...)
#		
#		end 
#
#	end 
#
#	ax.legend()
#
#end 
#
#
#
#
#
#
#
#
#

