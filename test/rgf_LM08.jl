include("rgf.jl")

function ribbon_x_armchair(NrAtomsPerLayer::Int)

	make_y_ribbon!(honeycomb_x_armchair(), NrAtomsPerLayer)

end 





fig,axes = PyPlot.subplots(1,3,num=8,figsize=(12,4))




L0 = make_y_ribbon!(Lattices.HoneycombLattice(),2)

for n=1:4

	plot_atoms(axes[1], Lattices.Atoms_ManyUCs(L0;Ns=n), label="n=$n")

	plot_atoms(axes[2], Lattices.Atoms_ManyUCs(ribbon_x_armchair(3),Ns=n),label="n=$n")

end

axes[1].legend() 
axes[2].legend() 

nr_layers = 16

for layer in 1:nr_layers

	atoms = slices_armchair(nr_layers, 12)[:AtomsOfLayer](layer)

	color = in(layer,1:nr_layers) ? (isodd(layer) ? "black" : "gray") : "red"

	plot_atoms(axes[3], atoms; color=color)

end 



















#ENERGIES = LinRange(0.001,0.2,50) 
#
#
#NrLayers =  360#N
#
#NrAtomsPerLayer = 1 + 3*23 # M 
#
#
#

















#
#
#
#
#
#
#
#
#function setind!(X,x,inds...)
#
#	@assert imag(x)<1e-10  x
#
#	X[inds...] = real(x)
#
#end 
#
#
#for (iE,Energy) in enumerate(ENERGIES)
#
#	gr = Gr(Energy)
##	ga = Ga(Energy) 
#
#	get_SE_ret = GreensFunctions.SelfEn_fromGDecim(gr, VirtLeads_ret, Slicer_ret)
##	get_SE_adv = GreensFunctions.SelfEn_fromGDecim(ga, VirtLeads_adv, Slicer_adv)
#
#	cc1 = ObservablesFromGF.CaroliConductance(gr, get_SE_ret, "Left","Right",2)
#
#	setind!(Y, cc1, 1, iE)
#
#	labels[1] = "T1"
#
#	cc2 = ObservablesFromGF.CaroliConductance2(gr, Ga(Energy), get_SE_ret, "Left","Right",1)
#
#	setind!(Y, cc2, 2, iE);labels[2] = "T2"
#
#
#
#	ff1 = ObservablesFromGF.FanoFactor2(gr, get_SE_ret, "Left", "Right", 2) 
#
#
#	setind!(Y,ff1,3,iE)
#
#	labels[3] = "FF"
#
##	Y[4,iE] = ObservablesFromGF.DOS(gr("Left",3,dir="left"))
#
#	labels[4] ="DOS"
#
#	print("\r",round(100iE/length(ENERGIES),digits=1),"%     ")
#
#
#end 
#
#println()#"\r                     ")
#
#PyPlot.figure(4)
#
#
#
#for (iy,(y,color)) in enumerate(zip(eachrow(Y),["blue","red","green","gold"]))
#
#
#	isapprox(y,zero(y),atol=1e-10) && continue
#
#	PyPlot.figure(4)
#
#	PyPlot.plot(ENERGIES,y,label=labels[iy],alpha=0.6,c=color) 
#
#	PyPlot.plot(-ENERGIES,y,alpha=0.6,c=color)
#
#
#
#end 
#
#PyPlot.legend()
#end 
