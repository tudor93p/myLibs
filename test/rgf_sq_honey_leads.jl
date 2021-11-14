include("rgf.jl")

PyPlot.close(101)
PyPlot.close(100)



nr_kpoints = 150 

dos_x = 0.6



fig1,Ax = PyPlot.subplots(2,3,num=100,figsize=(14,9),sharex=true,sharey=true) 

set_title(fig1)


	

fig2,Ax2 = PyPlot.subplots(1,1,num=101,figsize=(13,8),sharex=true,sharey=true) 

set_title(fig2)

lead_Energies = LinRange(-4.2,4.2, 150) 
	
lead_gf_labels = ["+","bulk"] 
	
DOS = zeros(2, 2, 2, length(lead_gf_labels), length(lead_Energies))

spectra = zeros(2, 2,2,nr_kpoints*4*nr_atoms_layer)

kLabels = LinRange(0,1,size(spectra)[end])

labels = fill("",2,2,2)


for (q,term_x) in enumerate([:armchair, :zigzag])

	LayerAtomRels = slices_ribbon(term_x, nr_layers, nr_atoms_layer)

	qAx = Ax[1:2, 2(q-1) .+ (term_x==:armchair ? [1,2] : [1])]


	nr_layers < 40 && plot_layers.(qAx; LayerAtomRels...)


	for i=1:2
		
		for j in (term_x==:armchair ? [1,2] : [1])
	
			ax = qAx[i,j] 
		
			for (args,k) in zip(two_lead_args(;LayerAtomRels...), (i,j))
		
				lab, get_latt = [(string(term_x), honeycomb_latt_lead),
												 ("square", square_latt_lead), 
																	 ][k]


				latt,a = get_latt(args...; LayerAtomRels...,four_uc=true,
												 term_x=term_x)
			
				label = string(args[3][1]," $lab")



				for n=0:1
		
					kw = n==0 ? Dict(:label=>label) : Dict()

					at = Lattices.Atoms_ManyUCs(latt; Ns=n)

					plot_atoms(ax, at; kw...,
										 c=only(colors[args[2].==[-1,1]])[n+1])
				end 
	
			end
		
			ax.legend()
	
				# q : armchair or zigzag 	
				#	i : square or honeycomb 
				#	j : left or right
	
	
			lab, get_latt = [(string(term_x), honeycomb_latt_lead),
												 ("square", square_latt_lead), 
																	 ][i]

			args = two_lead_args(;LayerAtomRels...)[j]

			term_x==:zigzag && lab=="square" && continue


			L1,a = get_latt(args...; LayerAtomRels...,four_uc=true, 
											term_x=term_x)
			
			label = string(args[3][1]," $lab")
	
	
	
	
			lead_hopp = get_lead_hopp(L1)
			
			labels[q,i,j] = label
	
	
			H = TBmodel.Bloch_Hamilt(Lattices.NearbyUCs(L1,1);
															 argH="k",
															 lead_hopp...)
		
			kPoints, kTicks = Utils.PathConnect(Lattices.BrillouinZone(L1), nr_kpoints; dim=2)
		
			spectrum = BandStructure.Diagonalize(H, kPoints; dim=2)
	
			spectra[q,i,j,:] = spectrum["Energy"]
	
	
			Ax2.set_xlabel(string("DOS / ",spectrum["kLabel"]," / DOS"))
	
			Ax2.set_ylabel("Energy")
		
			intra, inter = TBmodel.BlochHamilt_ParallPerp(L1, get_TBL, lead_hopp)
		
			intra_,inter_ = intra(), only(values(inter))
		
		
		#	i==1 && @show intra_ inter_ 
		
			for (ie,E) in enumerate(lead_Energies)
				
				gfs = GreensFunctions.GF_SanchoRubio(E+delta, intra_, inter_; 
																							 target=join(lead_gf_labels))
		
		#			@assert isapprox(GreensFunctions.GF_SanchoRubio(E + delta, intra_, inter_, "-"),
		#											 GreensFunctions.GF_SanchoRubio_rec(E + delta, intra_, inter_;																			 target="-")["-"],
		#											 rtol=1e-8)
		#
		#@assert isapprox(gfs["+"], l1[:GF](E)[1], rtol=1e-8)
		
		
				for (igf,gflab) in enumerate(lead_gf_labels)
		
					DOS[q,i,j,igf,ie] = ObservablesFromGF.DOS(gfs[gflab])
	
				end 
				
			
				print("\r$label ",(i,j)," ",round(100ie/length(lead_Energies),digits=1),"%     ") 
		
			end  #energy loop for lead of type i, on side j
		
			println()#"\r                     ")
	
	
		end #	j : left or right
	
	end #	i : square or honeycomb 

end # zigzag/armchair

sleep(0.5)

fig1.tight_layout()



#dev_TBi = ([0],[[0,0]],device_atoms(;LayerAtomRels...))
#
#dev_H = TBmodel.Bloch_Hamilt(dev_TBi; dev_hopp...)
#
#spectrum_device = BandStructure.Diagonalize(dev_H, zeros(0,1); dim=2)
#
#Ax2.scatter(spectrum_device["kLabels"], spectrum_device["Energy"], alpha=0.8, label= "Device",s=20,zorder=10,c="k") 






DOS = Utils.Rescale(DOS, [1,1+dos_x], [0, maximum(DOS)])  
	

shown_dos = []  

shown_spec = []


for q in 1:2, i in 1:2, j in 1:2 

	for (k,gflab) in enumerate(lead_gf_labels)

		d = DOS[q,i,j,k,:]

		isapprox(ones(size(d)),d,rtol=1e-7) && continue 

		any(shown_dos) do I
			
			isapprox(d,DOS[I...,:],rtol=1e-7)
			
		end && continue 

		label = string(labels[q,i,j]," ",lead_gf_labels[k])

		plot_d = k==1 ? d : -Utils.Rescale(d, [0,dos_x], [1,1+dos_x])

		 c=vcat(colors...)[length(shown_dos)+1]

		Ax2.plot(plot_d, lead_Energies; label=label,alpha=0.7, c=c)

		push!(shown_dos,(q,i,j,k))

	end  

	any(shown_spec) do I 

		isapprox(spectra[I...,:],spectra[q,i,j,:],rtol=1e-7)

	end && continue 


	isempty(labels[q,i,j]) || Ax2.scatter(kLabels,spectra[q,i,j,:], alpha=0.7,
							label=labels[q,i,j], s=10, zorder=length(shown_spec))

	push!(shown_spec, (q,i,j))

end 
	
Ax2.set_xlim(-dos_x,1+dos_x)
Ax2.legend() 
	




#fig1.tight_layout()
#fig2.tight_layout()

