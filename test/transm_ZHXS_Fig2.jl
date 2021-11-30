include("transm_ZHXS.jl")

f = ZHXS_hopping(0,0,0)[:Hopping]

for r in ([0,0],[1,0],[0,1],[2,1])

	break

	@show r
	
	println.(eachrow(f([0,0],r)))

	println()

end





PyPlot.close(2)

fig,Ax = PyPlot.subplots(2,2,num=2,figsize=(10,7))
sleep(0.1)

function figure2a(ax1,ax2,nr_at=100)

	nr_kpoints = 500

	hopping = ZHXS_hopping(0.35, -0.5, 1.5)

	ax2.set_title(string("N=",hopping[:NrMZM]))


	latt = device_lattice(nr_at)

	
	for n = -15:15

		plot_atoms(ax1, Lattices.Atoms_ManyUCs(latt; Ns=n), #label=n,
							 max_atoms=30,
							 color=["k","gray"][1+(n+100)%2],s=4)

	end 


	
	H = TBmodel.Bloch_Hamilt(Lattices.NearbyUCs(latt,1); 
													 argH="k", hopping...)


	y_oper = Operators.Position(2, Lattices.PosAtoms(latt); dim=2, hopping...)


	ky = Lattices.BrillouinZone(latt)[:,2]

	kPoints,kTicks = Utils.PathConnect(hcat(-ky/2,0*ky,ky/2), nr_kpoints; dim=2)


	spectrum = BandStructure.Diagonalize(H, kPoints; dim=2,
																			 kTicks=kTicks,
																			 operators = [["y"],[y_oper]])

	ax2.scatter(spectrum["kLabels"],spectrum["Energy"], s=5,
							c=spectrum["y"], cmap="viridis") 

	ax2.set_xlabel(spectrum["kLabel"])

	ax2.set_xticks(spectrum["kTicks"])
	ax2.set_xticklabels(["-π/2","0","π/2"])


	ax2.set_ylim(-.2,.2)

	ax2.set_xlim(0,1)

end 

function figure2bcd(ax,ax2,nr_at=100)

	latt = device_lattice(nr_at)

	NearbyUCs = Lattices.NearbyUCs(latt,1)
	ky = Lattices.BrillouinZone(latt)[:,2]
	
	kPoints,kTicks= Utils.PathConnect(hcat(-ky/2,0*ky,ky/2), 500; dim=2)


	atoms = Lattices.PosAtoms(latt)



	for (imu,mu) in enumerate([0, 0.6, 1.7])

		hopping = ZHXS_hopping(0.35, -0.5, mu) 
	
		ldos_oper = Operators.LDOS(; dim=2, nr_at=nr_at, nr_orb=hopping[:nr_orb])
		y_oper = Operators.Position(2, atoms; dim=2, nr_orb=hopping[:nr_orb])
		
		H = TBmodel.Bloch_Hamilt(NearbyUCs; argH="k", hopping...)
	
		spectrum = BandStructure.Diagonalize(H, kPoints; dim=2,
																				 operators=[["Y"],[y_oper]])

		spectrum0 = BandStructure.Diagonalize(H, zeros(2,1); dim=2,
																				 operators=[["LDOS"],[ldos_oper]])
	

		I = -0.2 .< spectrum["Energy"] .< 0.2 

		cmap = ["viridis","coolwarm","cool"][imu]

		ax2.scatter(spectrum["kLabels"][I], spectrum["Energy"][I],label ="$cmap:$mu", s=3, c=spectrum["Y"][I],cmap=cmap)

		ax2.set_xlabel(spectrum["kLabel"])

		N = hopping[:NrMZM]

		J = partialsortperm(abs.(spectrum0["Energy"]),1:(2N))

		ldos = sum(spectrum0["LDOS"][:,J], dims=2)

		ax.plot(atoms[2,:], sqrt.(ldos), label=string("mu=$mu, N=",N))

	end 

		ax2.set_xticks(kTicks)
		ax2.set_xticklabels(["-π/2","0","π/2"])

		ax2.legend()

	ax2.set_ylim(-.2,.2)
	ax.legend()

end 

	
figure2a(Ax[1,:]...,)

sleep(0.1) 

figure2bcd(Ax[2,:]...,)
	
	

