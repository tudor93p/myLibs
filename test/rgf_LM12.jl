include("rgf.jl")



PyPlot.close(121)
PyPlot.close(122)
PyPlot.close(123)


fig1,Ax1 = PyPlot.subplots(2,2,num=121,figsize=(10,9),sharex=true,sharey=true)


sleep(0.1)
fig3,Ax3 = PyPlot.subplots(2,2,num=123,figsize=(10,9),sharex=true,sharey=true)
sleep(0.1)
fig2,Ax2 = PyPlot.subplots(2,2,num=122,figsize=(8,7),sharex=true)

sleep(0.1) 

set_title(fig1)
set_title(fig2)





#plot_bonds(Ax[1,1], Rs_bonds; c="k",alpha=0.5)





svt_En = [chosen_energies["armchair"],chosen_energies["zigzag"]]


ENERGIES = sort(unique(vcat(svt_En...,
														Utils.uniqlinsp(-0.1,0.6,70, 3; Trunc=true))))

svt_En_ind = [[argmin(abs.(e .- ENERGIES)) for e in E] for E in svt_En]









for (i_latt_type,latt_type) in enumerate([:armchair,:zigzag])

	LayerAtomRels = slices_ribbon(latt_type, nr_layers, nr_atoms_layer)  

	inds_bonds,Rs_bonds = get_Bonds(;LayerAtomRels...)

	dev_at = device_atoms(;LayerAtomRels...)



	Labels,LeadLatts,LeadContacts,LHopps,CHopps = Utils.zipmap(
					two_lead_args(;LayerAtomRels...)) do (contact,dir,lab)
	
		latt,a = honeycomb_latt_lead(contact, dir, lab; 
															 LayerAtomRels..., term_x=latt_type,
#															 four_uc=true,
															 )
	
		lc = lead_contacts(contact, latt, a; LayerAtomRels...)

		lead_hopp = get_lead_hopp(latt)

		coupl_hopp = get_coupl_hopp(latt)
	
		return lab,latt,lc,lead_hopp,coupl_hopp
	
	end .|> collect
	
	
	leads_ret = prep_lead.(Labels, LeadLatts, LHopps, CHopps, delta)
	leads_adv = prep_lead.(Labels, LeadLatts, LHopps, CHopps, -delta)
	
	NG_ret = NewGeometry(LayerAtomRels, LeadContacts, leads_ret)
	
	G_ret = GF(NG_ret)  

	G_adv = GF(NewGeometry(LayerAtomRels, LeadContacts, leads_adv))
	

	y = zero(ENERGIES)

	y2 = zero(ENERGIES)





	lead_dos = zeros(2,3,length(ENERGIES))

	dos_labs = fill("",2,3)


#	SVT = zeros(2, 2, nr_layers*nr_atoms_layer)


	gfs = lead_gfs.(LeadLatts, [["bulk"]],delta)

	for (iE,Energy) in enumerate(ENERGIES)
	
		gr = G_ret(Energy)
		ga = G_adv(Energy)

		get_SE_ret = get_SE(gr, NG_ret) 
	
		cc1 = ObservablesFromGF.CaroliConductance(gr, get_SE_ret, Labels...)

		setIndRe!(y, cc1, iE)

		cc2 = ObservablesFromGF.CaroliConductance(gr,ga, get_SE_ret,
																							 reverse(Labels)...) 

		setIndRe!(y2, cc2, iE)


		for (i_lead,lab) in enumerate(Labels) 
			
			for uc in 1:2 

				lead_dos[i_lead,uc,iE] = ObservablesFromGF.DOS(gr(lab,uc,dir=lowercase(lab)))
			
				dos_labs[i_lead,uc] = "DOS $lab($uc)"

			end 

			uc= 3


			lead_dos[i_lead,uc,iE] = ObservablesFromGF.DOS(gfs[i_lead](Energy)[1])
			
			dos_labs[i_lead,uc] = "bulk $lab" 

		end 









		print("\r$latt_type leads, L=$nr_layers, W=$nr_atoms_layer: " ,
							round(100iE/length(ENERGIES),digits=1),"%     ")
		


		i_svt = indexin(iE,svt_En_ind[i_latt_type])[1]

		isnothing(i_svt) && continue 

		color=colors[i_svt][i_latt_type]

		Ax2[1,i_latt_type].scatter(Energy, y[iE], c=color,zorder=10,s=50)
	
		Ax2[2,i_latt_type].plot([Energy,Energy],[0,1],c=color, zorder=-10,linestyle=":",alpha=0.5)

		Ax2[1,i_latt_type].plot([Energy,Energy],[0,5],c=color, zorder=-10,linestyle=":",alpha=0.5)

		SVTs= [ObservablesFromGF.SiteTransmission(
												gr, ga,
												[ones(ComplexF64,1,1) for b in inds_bonds],
												inds_bonds,
												Rs_bonds,
												get_SE_ret,
												L,2; 
												dim=2) for L in Labels]


		#SVTs= [ObservablesFromGF.SiteTransmission0(
		#										gr, 
		#										[ones(ComplexF64,1,1) for b in inds_bonds],
		#										inds_bonds,
		#										Rs_bonds;
		#										dim=2) ]


		LDOS = ObservablesFromGF.LDOS_Decimation(gr; dim=2, LayerAtomRels...)


		ax = ax1,ax3 = [Ax[i_svt, i_latt_type] for Ax in [Ax1,Ax3]]
		

		for lead in leads_ret, a in ax 

				plot_atoms(a, lead[:head][1]; zorder=0, s=10,c="darkviolet",alpha=0.2)

		end



		plot_layers(ax1; LayerAtomRels..., c="k",zorder=0,s=5,
								max_nr_layers=nr_layers+1)



		for a in ax 
			a.set_title(string("$latt_type, E=",round(Energy,digits=3)))

		end 

#		for (iL,SVT) in enumerate((+(SVTs...),))
		for (iL,SVT) in enumerate(SVTs)

			color = colors[i_svt][i_latt_type+2(iL-1)]

			any(LinearAlgebra.norm.(eachcol(SVT)) .> 1e-5) || continue 

			ax1.quiver(eachrow(dev_at)..., eachrow(SVT)...;
							color=color, zorder=10)

		end  

		ax3.scatter(eachrow(dev_at)..., c=LDOS,cmap="PuBuGn",vmin=0,vmax=maximum(LDOS))

	end 
	
	println()




	Ax2[1,i_latt_type].plot(ENERGIES, y, c="k",alpha=0.5,zorder=-2)

#	isapprox(y,y2) 

	Ax2[1,i_latt_type].plot(ENERGIES, y2, c="red",alpha=0.5,zorder=-3)

	Ax2[1,i_latt_type].set_ylim(0,5) 

	Ax2[1,i_latt_type].set_xlim(extrema(ENERGIES))


	Ax2[1,i_latt_type].set_title("Conductance $latt_type")

	norm_dos = maximum(lead_dos[:,:,abs.(ENERGIES) .>= 0.05])

	lead_dos = reshape(lead_dos/norm_dos,:,length(ENERGIES))

	dos_labs = reshape(dos_labs,:)

	same = isapprox.(Algebra.OuterDist(lead_dos,lead_dos;dim=1)/length(ENERGIES),0,atol=1e-4)


	index_color = 1 
	for (r,dos) in enumerate(eachrow(lead_dos))

		any(same[1:r-1,r]) && continue 

		Ax2[2,i_latt_type].plot(ENERGIES, dos, label=dos_labs[r],
														c=colors[2][index_color],alpha=0.7)

		index_color +=1

	end 
	
	
	
	
	LS = lead_spectrum.(LeadLatts, 200)

	for (i_lead,spectrum) in enumerate(LS)

		E = spectrum["Energy"]

		any((isapprox(E, ls["Energy"]) for ls in LS[1:i_lead-1])) && continue 

		I = minimum(ENERGIES).<=E.<=maximum(ENERGIES)

		Ax2[2,i_latt_type].scatter(E[I], spectrum["kLabels"][I], label="PBC "*Labels[i_lead],s=10,c=colors[1][i_lead])

	end
		

	Ax2[2,i_latt_type].set_ylim(0,1)

	Ax2[2,i_latt_type].legend() 
end 









for fig in (fig1,fig2,fig3)

	fig.tight_layout()
	fig.subplots_adjust(top=0.92,hspace=0.1)

end 
