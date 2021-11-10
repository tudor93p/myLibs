include("rgf.jl")


PyPlot.close(11)







nr_layers_ = [60,120,180,300]

nr_atoms_layer_ = 3*[4, 8, 12, 24, 32, 36, 40, 49, 66].+1 

nr_layers0 = 180

#nr_layers_=nr_layers_[1:2]
#nr_atoms_layer_=nr_atoms_layer_[1:2]

ENERGIES = LinRange(0.001,2.7,100)

nr_obs = 2

Y = zeros(2, length(nr_layers_), length(nr_atoms_layer_),nr_obs,length(ENERGIES))

ylabs = fill("", nr_obs)


fig,Ax = PyPlot.subplots(3,nr_obs,num=11,figsize=(4*nr_obs,8),
												 #sharey="col"
												 )

for i_latt_type  in 1:2 

	latt_type, get_latt = [("square", square_latt_lead),
										 ("honey", honeycomb_latt_lead)
										 ][i_latt_type]


	for (i_len,Nr_Layers) in enumerate(nr_layers_)

		for (i_width,Nr_At_Layer) in enumerate(nr_atoms_layer_)
	
			layer_atom_rels = slices_armchair(Nr_Layers, Nr_At_Layer) 
		
			Labels,LeadLatts,LeadContacts,LHopps = Utils.zipmap(
							two_lead_args(;layer_atom_rels...)) do (contact,dir,lab)
		
				latt,a = get_latt(contact, dir, lab; layer_atom_rels...)
		
				lc = lead_contacts(contact, latt, a; layer_atom_rels...)
		
				lead_hopp = get_lead_hopp(latt)
		
				@assert length(lc)==Nr_At_Layer 
		
				return lab,latt,lc,lead_hopp
		
			end .|> collect
		
	
			leads_ret = prep_lead.(Labels, LeadLatts, LHopps, delta)
			leads_adv = prep_lead.(Labels, LeadLatts, LHopps, -delta)
		
			NG_ret = NewGeometry(layer_atom_rels, LeadContacts, leads_ret)
		
			G_ret = GF(NG_ret)  

			G_adv = GF(NewGeometry(layer_atom_rels, LeadContacts, leads_adv))
	
	
			for (iE,Energy) in enumerate(ENERGIES)
		
				latt_type!="square" && Nr_Layers!=nr_layers0 && break 

				gr = G_ret(Energy)
			
				get_SE_ret = get_SE(gr, NG_ret) 
			
				cc1 = ObservablesFromGF.CaroliConductance(gr, 
																									get_SE_ret, 
																									Labels..., 2)
		

				setIndRe!(Y, cc1*Nr_Layers/Nr_At_Layer, i_latt_type, i_len, i_width, 1, iE);  

				ylabs[1] = "Conductance"
	
	
	
				ff1 = ObservablesFromGF.FanoFactor2(gr, get_SE_ret, Labels..., 2)
			
				setIndRe!(Y,ff1,i_latt_type,i_len, i_width,2,iE)
				ylabs[2] = "FanoFactor"
			
			
				print("\r$latt_type leads, L=$Nr_Layers, W=$Nr_At_Layer: " ,
							round(100iE/length(ENERGIES),digits=1),"%     ")
		
				Nr_Layers!=nr_layers0 && break  

			end 
	
			println()
	
			for i_obs in 1:nr_obs

				if Nr_Layers==nr_layers0 

					Ax[i_latt_type,i_obs].plot(xy_PHS(ENERGIES, Y[i_latt_type,i_len,i_width,i_obs,:])...)  

				end 

				if latt_type=="square" 

					kw = i_width==1 ? Dict(:label=>Nr_Layers,) : ()

					Ax[3,i_obs].scatter(Nr_At_Layer/Nr_Layers,Y[1,i_len,i_width,i_obs,1];c=vcat(colors...)[i_len],s=5,kw...)

					Ax[3,i_obs].set_xlabel("W/L")
	
				end 
	
			end  
	
	
	
		end 
				


		for (i_obs,ylab) in enumerate(ylabs)

			Nr_Layers==nr_layers0 || continue 

			Ax[i_latt_type,i_obs].set_title("$latt_type leads, $ylab, L=$Nr_Layers")
		
			Ax[3,i_obs].set_title("square leads, $ylab")  
		
		end 
	end 

end 


[ax.legend() for ax in Ax[3,:]]


























fig.suptitle(string("Figure ",fig.number))


fig.tight_layout()
fig.subplots_adjust(top=0.95)
