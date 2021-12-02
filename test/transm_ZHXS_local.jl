include("transm_ZHXS.jl")
PyPlot.close(101) 
PyPlot.close(102) 


fig1,Ax1 = PyPlot.subplots(2,4,num=101,figsize=(20,12),sharex="col",sharey="col")

fig1.tight_layout()

sleep(0.1)



function figure101(ax, ax_)

	nr_layers = 20 

	widths = [80,150]

#	widths = div.(widths,2)

	widths[1] = 40

	BasisInfo = HoppingTerms.basis_info(HamiltBasis(true,true))
	dig3(x) = rstrip(rstrip(string(round(x,digits=3)),'0'),'.')


	mu_max = 1.7
	mu_min = 0 

	MUS0 = LinRange(mu_min,mu_max,5)[[2,3,4]]


	MUS0 = MUS0 .+ (rand(3).- 0.5)*(MUS0[2]-MUS0[1])*0.7



	#MUS0 = [0.15,0.66,1.6].+(2rand(3).-1)*0.001


	MUS0 = round.(MUS0,digits=3)

	MUS0 = [0.474, 0.813, 1.378]

	@show MUS0

	mu_min = min(mu_min,MUS0...)
	mu_max = max(mu_max,MUS0...)
	

	nr_MUS = 23 

	MUS = Utils.uniqlinsp(mu_min, mu_max, nr_MUS, 3; Trunc=true)


	for mu in vcat(MUS0...)

		MUS[argmin(abs.(MUS .- mu))] = mu
	
	end 
	


	iMUS0 = indexin(MUS0, MUS) 

	@assert all(!isnothing, iMUS0)

	delta = 1e-4

	Delta = 0.35  

	m = -0.5 

	i_w = 1 

	nr_at =40 

		println() 

		@show nr_at 



		ax2 = ax[2]

		ax4 = ax_[2]


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

		for a in [ax2,ax4]

			a.plot(fill(sqrt(m^2 -Delta^2),2),[0,1.1],color="green",ls="--",alpha=0.5)
			a.plot(extrema(MUS), fill(1/4,2),color="red",ls="--",alpha=0.5)

		end
	
	
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

				if length(MUS)<=50 
					
					for a in [ax2,ax4]
						a.scatter(mu,t[imu],color=col, alpha=0.8,s=4)
					end
				end 


				if imu>1 
	
					I = imu-1:imu  

					kw = Dict{Symbol,Any}(:color=>col,:alpha=>0.8)

					imu==2 && setindex!(kw, lab, :label)

					[a.plot(MUS[I], t[I]; kw...) for a in [ax2,ax4]]

				end  

			end 
	

			k = only(indexin(imu,iMUS0))
			
			if !isnothing(k)

				col = ["forestgreen","deeppink","olive"][k]

				ax1 = ax[[1,3,4][k]]
				ax3 = ax_[[1,3,4][k]]

				[a.plot([mu,mu],[0,1.1],lw=1,c=col,alpha=0.8,linestyle=":") for a in [ax2,ax4]]

				for t in (T,Tc,Tl)

					[a.scatter(mu,t[imu],c=col,s=30,zorder=10) for a in [ax2,ax4]]

				end
				
				
				alpha = 0.3

#				ax1.set_title(title)

				[LayeredSystem.plot_layers(a; LAR..., s=3,alpha=alpha,zorder=0) for a in [ax1,ax3]]

				for ll in lead_latts, n = 0:1, a in [ax1,ax3]

					Lattices.plot_atoms(Lattices.Atoms_ManyUCs(ll; Ns=n),a;s=3,max_atoms=51,alpha=alpha,color=["peru","dodgerblue"][isodd(n)+1],zorder=0)
	
				end  
			
			
				siteT = zeros(2,size(dev_at,2))



				for layer in 1:LAR[:NrLayers]


					IJ,XY,dXY = transversal_current(G_ret, G_adv, 
																			 GreensFunctions.DecayWidth∘se,
																			 layer, A; dev_hopping..., LAR..., dim=2)

					for ((i,j),t) in zip(eachcol(IJ),eachcol(dXY))

						siteT[:,i] += t 
						siteT[:,j] += t 

					end 


					layer in [1,LAR[:NrLayers]] && continue 


					dev_hopp_mat = Utils.add_args_kwargs(TBmodel.HoppingMatrix;
																						 dev_hopping...)


					IJ,XY,dXY = longitudinal_current(G_ret, G_adv, layer, dev_hopp_mat, NG_ret[2]; LAR..., dim=2)

					plot_arrows(ax1, XY, dXY; color=col, zorder=20, min_abs_len=0.1)



					if layer==div(LAR[:NrLayers],2) 
					
						s = sum(dXY)/2 

						println("\nNet longitudinal: 2*",dig3(s))
						
						for (t,lab,col) in [(T,"T","blue"),
															(Tc,"T_C","k"),(Tl,"T_L","red")] 

							println("$lab: ",dig3(t[imu]))

						end 

						plot_arrows(ax1, XY, dXY; color="k", zorder=30, min_abs_len=0.1)


						V = Utils.vectors_of_integers(3,1)
						u = [T[imu],Tc[imu],Tl[imu]]


						best_vs = V[partialsortperm(abs.(V*u .- s), 1:2),:]

						for v in eachrow(best_vs)

							q = findall(v.!=0) 

							isempty(q) && continue 

							str = join(string(v[qi] > 0 ? "+" : "-", ["T","TC","TL"][qi]) for qi in q)

							println("*) ",lstrip(str,'+'),"=",	dig3(LinearAlgebra.dot(v,u)))

						
						end 

						#println()




					end 
						

					for ((i,j),t) in zip(eachcol(IJ),eachcol(dXY))

						siteT[:,i] += t 
						siteT[:,j] += t 

					end 

#					plot_arrows(ax1, XY, dXY; color=col, zorder=20)

				end  



#				I = sort(Utils.flatmap(LAR[:IndsAtomsOfLayer], 2:LAR[:NrLayers]-1))


				i_mid_slice = LAR[:IndsAtomsOfLayer](div(LAR[:NrLayers],2)) 
				i_prev_mid_slice = LAR[:IndsAtomsOfLayer](div(LAR[:NrLayers],2)-1) 



				ax1.set_title("transm_ZHXS")

#				plot_arrows(ax1, dev_at, siteT; color=col, zorder=20,inds=I)

				ax3.set_title("ObservablesFromGF")

#				siteT3 = ObservablesFromGF.SiteTransmission(G_ret, hoppings_bonds, inds_bonds, Rs_bonds, se, A, 2; dim=2) 

#				plot_arrows(ax3, dev_at, siteT3; color=col, zorder=20,inds=I)


				only_horiz = [Ra[2]≈Rb[2] for (Ra,Rb) in Rs_bonds]

				bondT3 = -ObservablesFromGF.BondTransmission(G_ret, hoppings_bonds, inds_bonds, se, B, 2; dim=2) 


				inds_bonds_mid = map(inds_bonds) do (i,j)

					for (i0,j0) in zip(i_prev_mid_slice,i_mid_slice)

						i0==i && j0==j && return true 

					end 

					return false 

				end 







#				bondT3.*only_horiz.*[in(i,i_mid_slice)-in(j,i_mid_slice) for (i,j) in inds_bonds] |>sum |> println



				XY = mapreduce(Algebra.Mean, hcat, Rs_bonds)

				dXY = hcat((t*LinearAlgebra.normalize(Rb-Ra) for (t,(Ra,Rb)) in zip(bondT3,Rs_bonds))...)



				plot_arrows(ax3, XY, dXY; color=col, zorder=20, min_abs_len=0.1,inds=only_horiz)
				
				println("Net longitudinal method2: ",
								dig3(sum(dXY[1,inds_bonds_mid])))

				plot_arrows(ax3, XY, dXY; color="k", zorder=21, min_abs_len=0.1,inds=only_horiz.&inds_bonds_mid)
#				@show sum(siteT[1,i_mid_slice]) sum(siteT3[1,i_mid_slice])

				println()
			end 
	
			[a.set_ylim(0,1.1) for a in [ax2,ax4]]


			pad = 0.0

			[a.set_xlim(mu_min - pad*(mu_max-mu_min), mu_max + pad*(mu_max-mu_min)) for a in [ax2,ax4]]




			sleep(0.001)
			
			print("\r",Int(round(imu/length(MUS)*100)),"%           ")

		end 
		println()	
	
		[a.legend() for a in [ax2,ax4]]
	

#Ax2_[2,2].legend()

end  



figure101(Ax1[1,:],Ax1[2,:])

