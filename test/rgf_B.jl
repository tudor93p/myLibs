include("rgf.jl")


PyPlot.close(20)


fig1,Ax1 = PyPlot.subplots(1,1,num=20,figsize=(6,5))



set_title(fig1) 

sleep(0.05)


inds_bonds,Rs_bonds = get_Bonds(;LayerAtomRels..., NrLayers=10)



#plot_bonds(Ax1, Rs_bonds; c="k",alpha=0.5)

plot_layers(Ax1; LayerAtomRels..., c="k",zorder=0,s=10)

#dev_at = device_atoms(;LayerAtomRels...)




#q	= H_Superconductor.SC_Domain((ChemicalPotential=0,
#																Peierls=(0.1, 3sqrt(3)/2, :x),
#															 ),[1])[:Hopping] 
#
q = dev_hopp[:Hopping] 

cmap = "hsv" 

norm = PyPlot.cm.colors.Normalize(-1.0,1.0)   

col(Q::ComplexF64) = PyPlot.cm.get_cmap(cmap)(norm(angle(Q)/pi))



for Rs in Rs_bonds

	plot_bonds(Ax1, [Rs]; c=col(q(Rs...)[1]))

end 

fig1.colorbar(PyPlot.cm.ScalarMappable(norm=norm,cmap=cmap),ax=Ax1)























