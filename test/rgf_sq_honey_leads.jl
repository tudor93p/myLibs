include("rgf.jl")

PyPlot.close(101)
PyPlot.close(102)
PyPlot.close(100)



nr_kpoints = 150 

dos_x = 0.6



fig1,Ax = PyPlot.subplots(2,3,num=100,figsize=(12,8),sharex=true,sharey=true) 

set_title(fig1)

sleep(0.1)

fig2,Ax2 = PyPlot.subplots(1,1,num=101,figsize=(10,7),sharex=true,sharey=true) 


sleep(0.1)


set_title(fig2)

lead_Energies = LinRange(-4.2,4.2, 150) 
	
lead_gf_labels = ["+","bulk"] 
	
DOS = zeros(2, 2, 2, length(lead_gf_labels), length(lead_Energies))


						




four_uc = false#true

nr_eigval = Dict()

for (t,n) in [("square",four_uc ? 4 : 1),
							("zigzag", four_uc ? 8 : 4),
							("armchair", 4)]

	nr_eigval[t] = nr_kpoints*n*layer_width(Symbol(t), nr_atoms_layer)

end 


LDOS = zeros(2, 2, 2, maximum(length(get(chosen_energies,k,[])) for k in ("square","zigzag","armchair")), div(maximum(values(nr_eigval)),nr_kpoints))


fig3,Ax3 = PyPlot.subplots(3,size(LDOS,4),num=102,figsize=(12,4*size(LDOS,4)),sharex=true,sharey=true) 
sleep(0.1)



#LDOS[armchair/zigzag,sq/honey,left/right,ie,atoms]


LDOS_En_ind = Dict(k=>[argmin(abs.(lead_Energies .- e)) for e in get(chosen_energies,k,[])] for k in ("square","armchair","zigzag"))




spectra = zeros(2,2,2,maximum(values(nr_eigval)),2)


#LinRange(0,1,size(spectra)[end])

labels = fill("",2,2,2)


for (q,term_x) in enumerate([:armchair, :zigzag])

	LayerAtomRels = slices_ribbon(term_x, nr_layers, nr_atoms_layer)

	qAx = Ax[1:2, 2(q-1) .+ (term_x==:armchair ? [1,2] : [1])]


	nr_layers < 40 && plot_layers.(qAx; LayerAtomRels...)


	for i=1:2
			
		lab, get_latt = [(string(term_x), honeycomb_latt_lead),
												 ("square", square_latt_lead), 
																	 ][i]
		
		for j in (term_x==:armchair ? [1,2] : [1])

			ax = qAx[i,j] 
		
	
			for (args_,k_) in zip(two_lead_args(;LayerAtomRels...), (i,j))
			
				lab_, get_latt_ = [(string(term_x), honeycomb_latt_lead),
												 ("square", square_latt_lead), 
																	 ][k_]
	
	
				latt_,a_ = get_latt_(args_...; LayerAtomRels...,four_uc=four_uc,
												 term_x=term_x)
			
	
	
				for n_=0:1
			
					kw_ = n_==0 ? Dict(:label=>string(args_[3][1]," $lab_")) : ()

					at_ = Lattices.Atoms_ManyUCs(latt_; Ns=n_)
	
					plot_atoms(ax, at_; kw_...,
										 c=only(colors[args_[2].==[-1,1]])[n_+1])
				end 
		
			end
			
			ax.legend() 

				# q : armchair or zigzag 	
				#	i : square or honeycomb 
				#	j : left or right
	
		end 

		term_x==:zigzag && lab=="square" && continue


		for j=1:2

			args = two_lead_args(;LayerAtomRels...)[j]


			L1,a = get_latt(args...; LayerAtomRels...,four_uc=four_uc,
											term_x=term_x)
			
			label = string(args[3][1]," $lab")
	
	
	
			
			labels[q,i,j] = label
	
			spectrum = lead_spectrum(L1, nr_kpoints)	


			for (iks,ks) in enumerate(["kLabels","Energy"])
			
				spectra[q,i,j,1:nr_eigval[lab],iks] = spectrum[ks]

			end 
	
	
			Ax2.set_xlabel(string("DOS / ",spectrum["kLabel"]," / DOS"))
	
			Ax2.set_ylabel("Energy")

			gfs_ = lead_gfs(L1, lead_gf_labels, delta)

			for (ie,E) in enumerate(lead_Energies)

				gfs = gfs_(E)

				i_ldos = indexin(ie, LDOS_En_ind[lab])[1] 

				for (igf,gf) in enumerate(gfs)

					if !isnothing(i_ldos) && j!=igf 
						#&& igf==findfirst(isequal("bulk"),lead_gf_labels)

						nr_at = div(nr_eigval[lab],nr_kpoints)

						ldos = ObservablesFromGF.LDOS(gf;nr_at=nr_at,nr_orb=1)

						LDOS[q,i,j,i_ldos,1:nr_at] = ldos 

						DOS[q,i,j,igf,ie] = sum(ldos)
						
					else 

						DOS[q,i,j,igf,ie] = ObservablesFromGF.DOS(gf)

					end 

				end 
			

				print("\r$label ",(i,j)," ",round(100ie/length(lead_Energies),digits=1),"%     ") 
		
			end  #energy loop for lead of type i, on side j
		
			println()#"\r                     ")

			DOS[q,i,j,:,:] *= maximum(values(nr_eigval))/nr_eigval[lab]




		end #	j : left or right 



		for (igf,gflab) in enumerate(lead_gf_labels), j=1:2 

			igf==j && continue 

#			println((q,i,j)=>Dict("square"=>1,"armchair"=>2,"zigzag"=>3)[lab])

			args = two_lead_args(;LayerAtomRels...)[j]

			L1,a = get_latt(args...; LayerAtomRels..., four_uc=four_uc,
											term_x=term_x)

			for (i_ldos,e_ldos) in enumerate(chosen_energies[lab])
	
				ax3 = Ax3[Dict("square"=>1,"armchair"=>2,"zigzag"=>3)[lab],i_ldos]

	
				ax3.set_title(string("$lab, E=",round(e_ldos, digits=3),", bulk & $(args[2])"))

			
				ldos = LDOS[q,i,j,i_ldos,1:div(nr_eigval[lab],nr_kpoints)]#/maximum(LDOS[q,i,:,i_ldos,:])

				#println(i_ldos," ",maximum(ldos))

				for n in -div(nr_layers,8) .+2 .+ (gflab=="+" ? [0] : [0,1,2,3])
					
					plot_atoms(ax3, Lattices.Atoms_ManyUCs(L1; Ns=n);
										 c=ldos, cmap="PuBuGn",vmin=0,vmax=maximum(ldos),
										 s=7)

				end 

			end  #each energy 

		end # j:left or right 

		println()


	end #	i : square or honeycomb 

end # q: zigzag/armchair

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

		isapprox(spectra[I...,:,:], spectra[q,i,j,:,:], rtol=1e-7)

	end && continue 


	label = labels[q,i,j]

	if !isempty(label)

		N = only(v for (k,v) in nr_eigval if occursin(k,label))

		Ax2.scatter(
																				eachcol(spectra[q,i,j,1:N,:])...,
#																				kLabels,spectra[q,i,j,:], 
																				alpha=0.7,
						label=label, s=10, zorder=length(shown_spec))
	end 

	push!(shown_spec, (q,i,j))

end 
	
Ax2.set_xlim(-dos_x,1+dos_x)
Ax2.legend() 
	




fig3.tight_layout()
fig2.tight_layout() 
fig2.subplots_adjust(top=0.95)
#fig3.subplots_adjust(hspace=0.1)

