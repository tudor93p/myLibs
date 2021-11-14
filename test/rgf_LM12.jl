include("rgf.jl")



PyPlot.close(121)
PyPlot.close(122)


fig1,Ax1 = PyPlot.subplots(2,2,num=121,figsize=(9,8))
fig2,Ax2 = PyPlot.subplots(1,2,num=122,figsize=(9,5))

set_title(fig1)
set_title(fig2)

sleep(0.1)




#plot_bonds(Ax[1,1], Rs_bonds; c="k",alpha=0.5)





ENERGIES = LinRange(0.000,0.6,150)

svt_En = [0,0.3] 

svt_En_ind = [argmin(abs.(e .- ENERGIES)) for e in svt_En]









for (i_latt_type,latt_type) in enumerate([:armchair,:zigzag])


#	latt_type, get_latt = [("square", square_latt_lead),
#										 ("honey", honeycomb_latt_lead)
#										 ][3-i_latt_type]



	LayerAtomRels = slices_ribbon(latt_type, nr_layers, nr_atoms_layer)  

	inds_bonds,Rs_bonds = get_Bonds(;LayerAtomRels...)

	dev_at = device_atoms(;LayerAtomRels...)



	Labels,LeadLatts,LeadContacts,LHopps = Utils.zipmap(
					two_lead_args(;LayerAtomRels...)) do (contact,dir,lab)
	
	latt,a = honeycomb_latt_lead(contact, dir, lab; 
															 LayerAtomRels..., term_x=latt_type)
	
		lc = lead_contacts(contact, latt, a; LayerAtomRels...)
	
		lead_hopp = get_lead_hopp(latt)
	
		return lab,latt,lc,lead_hopp
	
	end .|> collect
	
	
	leads_ret = prep_lead.(Labels, LeadLatts, LHopps, delta)
	leads_adv = prep_lead.(Labels, LeadLatts, LHopps, -delta)
	
	NG_ret = NewGeometry(LayerAtomRels, LeadContacts, leads_ret)
	
	G_ret = GF(NG_ret)  

	G_adv = GF(NewGeometry(LayerAtomRels, LeadContacts, leads_adv))
	

	y = zero(ENERGIES)


#	SVT = zeros(2, 2, nr_layers*nr_atoms_layer)


	for (iE,Energy) in enumerate(ENERGIES)
	
		gr = G_ret(Energy)
	
		get_SE_ret = get_SE(gr, NG_ret) 
	
		cc1 = ObservablesFromGF.CaroliConductance(gr, get_SE_ret, Labels..., 2)

		setIndRe!(y, cc1, iE)


		print("\r$latt_type leads, L=$nr_layers, W=$nr_atoms_layer: " ,
							round(100iE/length(ENERGIES),digits=1),"%     ")
		


		i_svt = indexin(iE,svt_En_ind)[1]

		isnothing(i_svt) && continue 

		color=colors[i_svt][i_latt_type]

		Ax2[i_latt_type].scatter(Energy, y[iE], c=color,zorder=10,s=50)
		
		SVT= ObservablesFromGF.SiteTransmission(
												gr,
												[ones(ComplexF64,1,1) for b in inds_bonds],
												inds_bonds,
												Rs_bonds,
												get_SE_ret,
												Labels[1],2; 
												dim=2)


		ax = Ax1[i_svt, i_latt_type]
		
		
		ax.quiver(eachrow(dev_at)..., eachrow(SVT)...;
							color=color, zorder=10)

		plot_layers(ax; LayerAtomRels..., c="k",zorder=0,s=10)

		ax.set_title(string("$latt_type leads, E=",round(Energy,digits=2)))


	end 
	
	println()




	Ax2[i_latt_type].plot(ENERGIES, y, c="k",alpha=0.3,zorder=-2)

	Ax2[i_latt_type].set_ylim(0,5)
	Ax2[i_latt_type].set_xlim(extrema(ENERGIES))


	Ax2[i_latt_type].set_title("Conductance $latt_type leads")

end 












