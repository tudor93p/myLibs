include("transm_ZHXS.jl")
PyPlot.close(101) 
PyPlot.close(102) 


fig1,Ax1 = PyPlot.subplots(1,3,num=101,figsize=(14.5,6))

fig1.tight_layout()

sleep(0.1)



function figure101(ax, ax3=nothing)

	nr_layers = 20 

	widths = [80,150]

#	widths = div.(widths,2)

	widths[1] = 40

	BasisInfo = HoppingTerms.basis_info(HamiltBasis(true,true))

	MUS0 = [[0.15,1.6] for i=1:2]

	mu_min = min(0,minimum(minimum, MUS0))
	mu_max = max(1.7,maximum.(MUS0)...)
	

	nr_MUS = 20 

#	MUS = sort(unique(vcat(MUS0...,
#												 Utils.uniqlinsp(mu_min, mu_max, nr_MUS, 2; Trunc=true)
#												 )))

	MUS = Utils.uniqlinsp(mu_min, mu_max, nr_MUS, 2; Trunc=true)

	for mu in vcat(MUS0...)

		MUS[argmin(abs.(MUS .- mu))] = mu
	
	end 
	



	iMUS0 = indexin.(MUS0,[MUS]) 

	delta = 1e-4

	Delta = 0.35  

	m = -0.5 

	for (i_w, nr_at) in enumerate(widths) 

		i_w==1||continue 

		println() 

		@show nr_at 


#		ax1 = ax[i_w,[1,3]]

		ax2 = [ax,ax3][i_w][2]




		LAR = slices_ribbon(nr_layers, nr_at)  

		dev_at = device_atoms(;LAR...)

		inds_bonds, Rs_bonds = get_Bonds(;LAR...) 

		dev_hopping0 = ZHXS_hopping(Delta, m, 100rand()) 

		hoppings_bonds = [dev_hopping0[:Hopping](Rs...) for Rs in Rs_bonds]

#		@show rand(hoppings_bonds)
#		Lattices.plot_bonds(Rs_bonds, ax1; c="gray",lw=0.4)




	
		lead_latts,lead_labels,lead_conts = Utils.zipmap(two_lead_args(;LAR...)) do args 
	
			ll = square_latt_lead(args...; LAR...)
	
	
	
			return ll,args[3], lead_contacts(args[1], ll;  LAR...)
		
		end  .|> collect
	
	
	

	

#		title = string("(",'a'+i_w-1,") $nr_layers x $nr_at")
		title = string("$nr_layers x $nr_at")


#		ax2.set_title("$title, m=$m")
			

		T,Tl,Tc = [zeros(length(MUS)) for i=1:3]
	
		A, B = lead_labels 
	
		ax2.plot(fill(sqrt(m^2 -Delta^2),2),[0,1.1],color="green",ls="--",alpha=0.5)
		ax2.plot(extrema(MUS), fill(1/4,2),color="red",ls="--",alpha=0.5)
	
	
		lead_hopping = ZHXS_hopping(0, m, 0) 
		
		
		leads_ret = prep_lead.(lead_labels, lead_latts, [lead_hopping], delta) 
		leads_adv = prep_lead.(lead_labels, lead_latts, [lead_hopping], -delta)  
		



		el_proj, hole_proj = [Operators.Projector(BasisInfo(x), :atoms; nr_orb=length(BasisInfo(x)), nr_at=nr_at, dim=2) for x in "EH"]
	
		for (imu,mu) in enumerate(MUS)
	
			dev_hopping = ZHXS_hopping(Delta, m, mu) 
		
			NG_ret = NewGeometry(LAR, lead_conts, leads_ret, dev_hopping)
			G_ret = GF(dev_hopping, NG_ret)(0)
			se = get_SE(G_ret, NG_ret)

			NG_adv = NewGeometry(LAR, lead_conts, leads_adv, dev_hopping)
			
			G_adv = GF(dev_hopping, NG_adv)(0)
		
#			@assert G_ret(A,1,B,1)' ≈ G_adv(B,1,A,1)
		
	
		
	
			T[imu] = ObservablesFromGF.CaroliConductance(G_ret, se, A, B, el_proj)
	
			Tl[imu] = ObservablesFromGF.CaroliConductance(G_ret, se, A, el_proj, hole_proj)
	
	
			Tc[imu] = real(ObservablesFromGF.CaroliConductance(G_ret, G_adv, se, A, B, el_proj, hole_proj))

			for (t,lab,col) in [(T,"T","blue"),
															(Tc,"T_C","k"),(Tl,"T_L","red")]

				length(MUS)<=50&& ax2.scatter(mu,t[imu],color=col, alpha=0.8,s=4)

				if imu>1 
	
					I = imu-1:imu  

					kw = Dict{Symbol,Any}(:color=>col,:alpha=>0.8)

					imu==2 && setindex!(kw, lab, :label)

					ax2.plot(MUS[I], t[I]; kw...)

				end  

			end 
	

			k = only(indexin(imu,iMUS0[i_w]))
			
			if !isnothing(k)

				col = ["forestgreen","deeppink"][k]

				ax1 = [ax,ax3][i_w][[1,3][k]]

				ax2.plot([mu,mu],[0,1.1],lw=1,c=col,alpha=0.8,linestyle=":")

				ax2.scatter(mu,T[imu],c=col,s=10,zorder=10)
				
				
				alpha = 0.3

#				ax1.set_title(title)

				LayeredSystem.plot_layers(ax1; LAR..., s=3,alpha=alpha,zorder=0)

				for ll in lead_latts, n = 0:5

					Lattices.plot_atoms(Lattices.Atoms_ManyUCs(ll; Ns=n),ax1;s=3,max_atoms=51,alpha=alpha,color=["peru","dodgerblue"][isodd(n)+1],zorder=0)
	
				end  

				for layer in 1:LAR[:NrLayers]

					XY,dXY = transversal_current(G_ret, G_adv, 
																			 GreensFunctions.DecayWidth∘se,
																			 layer, A; dev_hopping..., LAR..., dim=2)
			
					for (xy,dxy) in zip(eachcol(XY),eachcol(dXY))

						plot_arrow(ax1, xy, dxy; color=col, zorder=20)

					end 

				end 

#				siteT = ObservablesFromGF.SiteTransmission(G_ret, hoppings_bonds, inds_bonds, Rs_bonds, se, A, 1; dim=2)
#
#				max_norm = maximum(LinearAlgebra.norm, eachcol(siteT))
#
#
#			for (xy,dxy) in zip(eachcol(dev_at),eachcol(siteT))
#
#				LinearAlgebra.norm(dxy)<0.07max_norm && continue 
#
#				#ax1.quiver(eachrow(dev_at)...,  eachrow(siteT)..., color=col,zorder=20,pivot="mid", length_includes_head=true, width=0.1, head_width=0.41, head_length=0.2)
#
#				plot_arrow(ax1, xy, dxy; color=col, zorder=20)
#
#			end


			end 
	
			ax2.set_ylim(0,1.1) 


			pad = 0.0

			ax2.set_xlim(mu_min - pad*(mu_max-mu_min), mu_max + pad*(mu_max-mu_min))




			sleep(0.001)
			
			print("\r",Int(round(imu/length(MUS)*100)),"%           ")

		end 
		println()	
	
		ax2.legend()
	end	
	

#Ax2_[2,2].legend()

end  



figure101(Ax1)

