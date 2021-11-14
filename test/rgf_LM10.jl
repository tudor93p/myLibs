include("rgf.jl")


PyPlot.close(10)


fig,Ax = PyPlot.subplots(2,2,num=10,figsize=(12,10),
												 sharex="col",sharey="col")

set_title(fig)




for i=1:2 

	ax = Ax[i,1]

	nr_layers<100 && plot_layers(ax; LayerAtomRels...) 

	lab, get_latt = [("square", square_latt_lead),
										 ("honeycomb", honeycomb_latt_lead)
															 ][i]


	Labels,LeadLatts,LeadContacts,LHopps = Utils.zipmap(
										[(1, -1, "Left",["orange","green"]), 
										(LayerAtomRels[:NrLayers], 1, "Right",["blue","gold"])]
								 ) do args

		latt,a = get_latt(args[1:3]...; LayerAtomRels...)

		label = string(args[3]," $lab")


		for n=0:1

			nr_layers<100 || continue

			kw = n==0 ? Dict(:label=>label) : Dict() 

			plot_atoms(ax, Lattices.Atoms_ManyUCs(latt; Ns=n); kw...,
								 color=args[4][n+1])
		
		end 

		lc = lead_contacts(args[1], latt, a; LayerAtomRels...)

		lead_hopp = get_lead_hopp(latt)

		@assert length(lc)==nr_atoms_layer 

	test1 = isapprox(1,lead_hopp[:Hopping]([sqrt(3),0],[0,0])[1]) && isapprox(1,lead_hopp[:Hopping]([0,0],[0,sqrt(3)])[1]) 

	test2 = isapprox(1,lead_hopp[:Hopping]([1,0],[0,0])[1]) && isapprox(1,lead_hopp[:Hopping]([0,0],[0,1])[1])

	@assert [test1,test2][i] [test1,test2]

		return label,latt,lc,lead_hopp

	end .|> collect

	ax.legend()


	leads_ret = prep_lead.(Labels, LeadLatts, LHopps, delta)
	leads_adv = prep_lead.(Labels, LeadLatts, LHopps, -delta)

	NG_ret = NewGeometry(LayerAtomRels, LeadContacts, leads_ret)
	NG_adv = NewGeometry(LayerAtomRels, LeadContacts, leads_adv)

	G_ret = GF(NG_ret) 
	G_adv = GF(NG_adv) 

	
	ax = Ax[i,2] 



	ENERGIES = LinRange(0.001,0.2,50) 

	Y = zeros(4,length(ENERGIES))

	ylabs = ["" for y in eachrow(Y)]

	for (iE,Energy) in enumerate(ENERGIES)
	
		gr = G_ret(Energy)
	
		get_SE_ret = get_SE(gr, NG_ret) 
	
		cc1 = ObservablesFromGF.CaroliConductance(gr, 
																							get_SE_ret, 
																							Labels..., 2)
	
		setIndRe!(Y, cc1, 1, iE); ylabs[1] = "T1"





		ga = G_adv(Energy)
		get_SE_adv = get_SE(ga, NG_adv) 

		A,B = Labels 
		uc=1 

		a1,a2 = rand(1:100,2)

		@assert isapprox(gr("Atom",a1,"Atom",a2),ga("Atom",a2,"Atom",a1)')


#		println(isapprox(ga(A,uc,B,uc),gr(B,uc,A,uc)',rtol=1e-6))

		cc2 = ObservablesFromGF.CaroliConductance2(gr, ga, get_SE_ret,
																							Labels..., 2)

		if !isapprox(cc1,cc2)

			setIndRe!(Y, cc2, 2, iE); ylabs[2] = "T2"

		end
	
	
		ff1 = ObservablesFromGF.FanoFactor2(gr, get_SE_ret, Labels..., 2)
	
	
		setIndRe!(Y,ff1,3,iE); ylabs[3] = "FF"
	
	#	Y[4,iE] = ObservablesFromGF.DOS(gr("Left",3,dir="left"))
	
#		labels[4] ="DOS"
	
		print("\r$i ",round(100iE/length(ENERGIES),digits=1),"%     ")
	
	end 

	println()

	nonzeros = map(zip(eachrow(Y),ylabs)) do (y,lab)
		
		isapprox(y,zero(y),atol=1e-10) && return false 
	
		return true 
		#!isempty(lab) || 

	end


	for (y,lab,c) in zip(eachrow(Y[nonzeros,:]),ylabs[nonzeros],
											 vcat(colors...)[1:count(nonzeros)])

		ax.plot(ENERGIES,y,label=lab,alpha=0.6,c=c) 
		ax.plot(-ENERGIES,y,alpha=0.6,c=c) 
	
	end 	
	
	ax.legend()

	ax.set_xlim(maximum(ENERGIES)*[-1,1])

#	ax.set_ylim(-3,5)#extrema(Y))


end 










