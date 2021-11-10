include("rgf.jl")


PyPlot.close(20)


fig1,Ax1 = PyPlot.subplots(1,1,num=20,figsize=(6,5))



set_title(fig1) 

sleep(0.05)


inds_bonds,Rs_bonds = get_Bonds(;LayerAtomRels...)



plot_bonds(Ax1, Rs_bonds; c="k",alpha=0.5)

plot_layers(Ax1; LayerAtomRels..., c="k",zorder=0,s=10)

#dev_at = device_atoms(;LayerAtomRels...)




