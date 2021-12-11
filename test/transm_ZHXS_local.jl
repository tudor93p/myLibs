include("transm_ZHXS.jl") 




function figure101(ax, ax_)

	nr_layers = 20 
	nr_at = rand(20:32)

	BasisInfo = HoppingTerms.basis_info(HamiltBasis(true,true))
	dig3(x) = rstrip(rstrip(string(round(x,digits=3)),'0'),'.')

	delta = 1e-4

	Delta = 0.35  

	m = -0.5 

	show_transversal = false 
	show_longitudinal = false
	show_site = true


	min_abs_len = 0.2

	mu_max = 1.7
	mu_min = 0

	MUS0 = LinRange(mu_min,mu_max,5)[[2,3,4]] 

	MUS0 = MUS0 .+ (rand(3).- 0.5)*(MUS0[2]-MUS0[1])*0.7


#	MUS0 = [0.15,0.66,1.6].+(2rand(3).-1)*0.01
	MUS0 = [sqrt(m^2 -Delta^2)/2.1, 1.356, 1.583]

	MUS0 .+= (2rand(length(MUS0)).-1)*0.01



	MUS0 = round.(MUS0, digits=3)


	@show MUS0

	mu_min = min(mu_min,MUS0...)
	mu_max = max(mu_max,MUS0...)
	

	nr_MUS = 10

	MUS = sort(unique(vcat(MUS0, Utils.uniqlinsp(mu_min, mu_max, nr_MUS, 3; Trunc=true))))


	for mu in MUS0

		MUS[argmin(abs.(MUS .- mu))] = mu
	
	end 
	


	iMUS0 = indexin(MUS0, MUS) 

	@assert all(!isnothing, iMUS0)



	println() 

#	@show nr_at 



	ax2 = ax[2]

	ax4 = ax_[2]


	LAR = slices_ribbon(nr_layers, nr_at)  

	dev_at = device_atoms(;LAR...)

	inds_bonds, Rs_bonds = get_Bonds(;LAR...) 

	midpoints_bonds = mapreduce(Algebra.Mean, hcat, Rs_bonds)


	slice_bonds_ = get_local_bonds(show_transversal, show_longitudinal,inds_bonds,Rs_bonds; LAR...)

	dev_hopping0 = ZHXS_hopping(Delta, m, 100rand()) 

	hoppings_bonds = [dev_hopping0[:Hopping](Rs...)' for Rs in Rs_bonds]

#	@show rand(hoppings_bonds)
#	Lattices.plot_bonds(Rs_bonds, ax1; c="gray",lw=0.4)




	
	lead_latts,lead_labels,lead_conts = Utils.zipmap(two_lead_args(;LAR...)) do args 
	
		ll = square_latt_lead(args...; LAR...)
	
		return ll,args[3], lead_contacts(args[1], ll;  LAR...)
	
	end  .|> collect
	
	
	

	

#	title = string("(",'a'+i_w-1,") $nr_layers x $nr_at")
#	title = string("$nr_layers x $nr_at")



#	ax2.set_title("$title, m=$m")
		

	T,Tl,Tc = [zeros(length(MUS)) for i=1:3]
	
	A, B = lead_labels 

	for a in [ax2]

		a.plot(fill(sqrt(m^2 -Delta^2),2),[0,1.1],color="green",ls="--",alpha=0.5,lw=0.5)
		a.plot(extrema(MUS), fill(1/4,2),color="red",ls="--",alpha=0.5,lw=0.5)

	end
	
	
	lead_hopping = ZHXS_hopping(0, m, 0) 
	
	
	leads_ret = prep_lead.(lead_labels, lead_latts, [lead_hopping], delta) 
	leads_adv = prep_lead.(lead_labels, lead_latts, [lead_hopping], -delta)  
#
#	leads1 = prep_lead.(lead_labels, lead_latts, [lead_hopping])
#	leads2 = prep_lead.(lead_labels, lead_latts, [lead_hopping])

#	Energy = rand()*0.001 + delta*im
#
#	gf_lead_ret = leads1[1][:GF](Energy)[1] 
#	gf_lead_adv = leads1[1][:GF](Energy')[1] 
#
#	@assert gf_lead_ret ≈ leads_ret[1][:GF](real(Energy))[1] 
#	@assert gf_lead_adv ≈ leads_adv[1][:GF](real(Energy))[1] 

	Energy = 0.0

	el_proj, hole_proj = [Operators.Projector(BasisInfo(x), :atoms; nr_orb=length(BasisInfo(x)), nr_at=nr_at, dim=2) for x in "EH"] 

#	projs_ = Dict("identity"=>identity,"hole"=>hole_proj,"electron"=>el_proj)

#	projs_k = collect(Base.product((sort(collect(keys(projs_))) for _=1:2)...))[i_projs]

#	@info "$i_projs) Projectors: "*join(projs_k,", ")

#	projs = [projs_[p] for p in projs_k]


	for (imu,mu) in enumerate(MUS)
	
		print("\r",Int(round((imu-1)/length(MUS)*100)),"%  mu=$mu         ")

		dev_hopping = ZHXS_hopping(Delta, m, mu) 

#		NG1 = NewGeometry(LAR, lead_conts, leads1, dev_hopping)  
#		NG2 = NewGeometry(LAR, lead_conts, leads2, dev_hopping)  
#
#		G_ret2 = GF(dev_hopping, NG1; leads_have_imag=false)(Energy)
#		G_adv2 = GF(dev_hopping, NG2; leads_have_imag=false)(Energy')

		NG_ret = NewGeometry(LAR, lead_conts, leads_ret, dev_hopping) 
		NG_adv = NewGeometry(LAR, lead_conts, leads_adv, dev_hopping)

		G_adv = GF(dev_hopping, NG_adv; leads_have_imag=true)(real(Energy)) 
		G_ret = GF(dev_hopping, NG_ret; leads_have_imag=true)(real(Energy)) 


#		@assert G_ret1(A,1,dir="left") ≈ gf_lead_ret
#		@assert G_adv1(A,1,dir="left") ≈ gf_lead_adv 
#
#		@assert G_adv2(A,1,dir="left") ≈ gf_lead_adv 
#		@assert G_ret2(A,1,dir="left") ≈ gf_lead_ret
#
#		@assert G_ret1(A,1,B,1)' ≈ G_adv1(B,1,A,1)
#		@assert G_ret2(A,1,B,1)' ≈ G_adv2(B,1,A,1)
#
#
#		se1,se2 = get_SE(G_ret1, NG_ret)(A), get_SE(G_ret2, NG1)(A) 
#
#		@assert se1≈se2 LinearAlgebra.norm(se1-se2)
#

		se = get_SE(G_ret, NG_ret) 



	
	
	
	
	
		T[imu] = ObservablesFromGF.CaroliConductance(G_ret, se, A, B, el_proj)
			
		Tl[imu] = ObservablesFromGF.CaroliConductance(G_ret, se, A, el_proj, hole_proj)
			
		Tc[imu] = ObservablesFromGF.CaroliConductance(G_ret, se, A, B, el_proj, hole_proj)
		
		@assert Tc[imu] ≈ real(ObservablesFromGF.CaroliConductance(G_ret, G_adv, se, A, B, el_proj, hole_proj))

#		ObservablesFromGF.DOS_Decimation(G_ret; LAR..., proj=el_proj, dim=2, NG_ret[4]...)|>println

#ObservablesFromGF.ComputeDOSLDOS_Decimation(G_ret, nr_layers, LAR[:IndsAtomsOfLayer], el_proj,NG_ret[4]; dos=true,ldos=true,dim=2) |> typeof |> println

#		continue 

		for (t,lab,col) in [(T,"T","blue"),
														(Tc,"T_C","k"),(Tl,"T_L","red")]

			if length(MUS)<=50 
				
				for a in [ax2]
					a.scatter(mu,t[imu],color=col, alpha=0.8,s=4)
				end
			end 


			if imu>1 
	
				I = imu-1:imu  

				kw = Dict{Symbol,Any}(:color=>col,:alpha=>0.8)

				imu==2 && setindex!(kw, lab, :label)

				[a.plot(MUS[I], t[I]; kw...) for a in [ax2]]

			end  

		end 


		k = only(indexin(imu,iMUS0))
		
		if !isnothing(k)

			col = ["forestgreen","deeppink","olive"][k]



			ax1 = ax[[1,3,4][k]]
			ax3 = ax_[[1,3,4][k]]

			[a.set_title("mu = $mu") for a in [ax1,ax3]] 

			[a.plot([mu,mu],[0,1.1],lw=2,c=col,alpha=0.8,linestyle=":") for a in [ax2]]

			for (t,lab,col_) in [(T,"T","blue"),
														(Tc,"T_C","k"),(Tl,"T_L","red")]

				ax2.scatter(mu,t[imu],c=col_,s=40,zorder=10)

				if show_longitudinal  

					ax4.scatter([2,3,nr_layers+1,nr_layers], fill(t[imu],4), c=col_,s=40,zorder=30, alpha=0.7)

				end 
					
				if show_site 
					
					ax4.scatter([0,1,nr_layers+1,nr_layers], fill(t[imu],4), c=col_,s=40,zorder=30, alpha=0.7)

				end 
			end
			
			
			alpha = 0.3

			[LayeredSystem.plot_layers(a; LAR..., s=3,alpha=alpha,zorder=0) for a in [ax1,ax3]]

			for ll in lead_latts, n = 0:1, a in [ax1,ax3]

				Lattices.plot_atoms(Lattices.Atoms_ManyUCs(ll; Ns=n),a;s=3,max_atoms=51,alpha=alpha,color=["peru","dodgerblue"][isodd(n)+1],zorder=0)
	
			end  
		
		
			siteT = zeros(2,size(dev_at,2))

			T3_lead = (B,1)
#			T3_lead = (A,1)
			

			bondT3 = ObservablesFromGF.BondTransmission(G_ret, G_adv, 
																									hoppings_bonds, 
																									inds_bonds, se, 
																									T3_lead...) 
		

			siteT3 = ObservablesFromGF.SiteTransmission(G_ret, hoppings_bonds, inds_bonds, Rs_bonds, se, T3_lead...; dim=2)

			@assert siteT3 ≈ ObservablesFromGF.SiteTransmission(G_ret, G_adv, hoppings_bonds, inds_bonds, Rs_bonds, se, T3_lead...; dim=2)

			siteT3 = ObservablesFromGF.SiteTransmission0(G_ret, hoppings_bonds, inds_bonds, Rs_bonds; dim=2)

			
			@assert bondT3 ≈ ObservablesFromGF.BondTransmission(G_ret, hoppings_bonds, inds_bonds, se, T3_lead...)

			bondT3_arrows = hcat((t*LinearAlgebra.normalize(Rb-Ra) for (t,(Ra,Rb)) in zip(bondT3,Rs_bonds))...)



			for layer in 1:nr_layers

				layer==1 && println()

				IJ,XY,dXY = transversal_current(G_ret, G_adv, 
																		 GreensFunctions.DecayWidth∘se,
																		 layer, T3_lead...; dev_hopping..., LAR..., dim=2)

				for ((i,j),t) in zip(eachcol(IJ),eachcol(dXY))

					siteT[:,i] += t 
					siteT[:,j] += t 

				end 

				if show_transversal  

					slice_bonds = slice_bonds_(layer) 
					
					plot_arrows(ax1, XY, dXY; color=col, zorder=20, min_abs_len=min_abs_len)

					XY2 = midpoints_bonds[:,slice_bonds]

					dXY2 = bondT3_arrows[:,slice_bonds] 

					plot_arrows(ax3, XY2, dXY2; color=col, zorder=20, min_abs_len=min_abs_len)


					@assert XY ≈ XY2 

					d = Algebra.Mean(abs, dXY2[1,:]-dXY[1,:])  

					text = "Layer $layer: $d" 

					
#					d<1e-3 && println(text)
					d<1e-3 || @warn text 
						

#					layer == 3 && error()

#					d > 1e-3 && @warn string("Mean transersal difference: ",dig3(d))
				end

				1<layer<nr_layers || continue  




				slice_bonds = slice_bonds_(layer) 



				dev_hopp_mat = Utils.add_args_kwargs(TBmodel.HoppingMatrix;
																					 dev_hopping...)





				IJ,XY,dXY = longitudinal_current(G_ret, G_adv, layer, dev_hopp_mat, NG_ret[2]; LAR..., dim=2) 


				if show_longitudinal 
					
					plot_arrows(ax1, XY, dXY; color=col, zorder=20, min_abs_len=min_abs_len)
					plot_arrows(ax3, midpoints_bonds, bondT3_arrows; color=col, zorder=20, min_abs_len=min_abs_len, inds=slice_bonds) 


					d = Algebra.Mean(abs, bondT3_arrows[1,slice_bonds]-dXY[1,:]) 
#					d > 1e-3 && @warn string("Mean longitudinal difference: ",dig3(d))

				end



				for ((i,j),t) in zip(eachcol(IJ),eachcol(dXY))

					siteT[:,i] += t 
					siteT[:,j] += t 

				end 




				s = sum(dXY[1,:])/2 


				if show_longitudinal 

#					println("  $layer) Net longitudinal method 1: 2*",dig3(s)) 


					ax4.scatter(layer,s,c=col,zorder=10,alpha=0.7,s=30)


				end 
				
				if layer==div(nr_layers,2) 

					for (t,lab,col) in [(T,"T","blue"),
														(Tc,"T_C","k"),(Tl,"T_L","red")] 
	
						println("$lab: ",dig3(t[imu]))
	
					end 

					X = dXY[1,:]

					if show_longitudinal 
						plot_arrows(ax1, XY, dXY; color="k", zorder=30, min_abs_len=min_abs_len)
	
#						println("Black 1: [",join(dig3.(X[abs.(X).>0.1]),", "),"]")


						for (trial,v) in zip(best_linear_comb(s, 1, 2, T[imu], Tc[imu], Tl[imu])...)
							q = findall(v.!=0) 
		
							str = join(string(v[qi] > 0 ? "+" : "-", ["T","TC","TL"][qi]) for qi in q)
	
							println("*) ",lstrip(str,'+'),"=",	dig3(trial))
							
						end 
					end 
	

				end # if layer == middle 









				
				
				if show_longitudinal 

					X = bondT3_arrows[1,slice_bonds]
				
	

					s2 = sum(X)/2
	
					pr_s2 = string("  $layer) Net longitudinal method 2: 2*", dig3(s2)) 

#					println(pr_s2)

					#abs(s-s2)/maximum(abs,(s,s2)) > 1e-2 && @warn pr_s2 
					
	
					if layer == div(nr_layers,2)
	
	
#						println("Black 2: [",join(dig3.(X[abs.(X).>0.1]),", "),"]")
	
						plot_arrows(ax3, midpoints_bonds, bondT3_arrows; 
												color="k", zorder=21, min_abs_len=min_abs_len, inds=slice_bonds)
	
					end 

				end 


			end   # for each layer 

			if show_longitudinal	
				x = 2:nr_layers 
	
				ax4.plot(x,
								 [sum(bondT3_arrows[1,slice_bonds_(L)])/2 for L in x],
								 color=col,zorder=5,alpha=0.7)
	
				ax4.set_xlim(1,nr_layers+1)

			end 

			ax2.set_title("transm_ZHXS")


			ax4.set_title("ObservablesFromGF")

			if show_site 
				
				J = mapreduce(LAR[:IndsAtomsOfLayer], vcat, 2:nr_layers)

				println("Mean site diff: ",
								dig3(Algebra.Mean(LinearAlgebra.norm, 
																	collect(eachcol(siteT-siteT3))[J])))





				for layer in 1:nr_layers 
	
					I = LAR[:IndsAtomsOfLayer](layer)   

					if 2<layer<nr_layers-1

						sf = 1/4
	
						ax4.scatter(layer,sum(siteT[1,I])*sf,color=col,zorder=5,alpha=0.7)
				
						s1 = sum(siteT3[1, LAR[:IndsAtomsOfLayer](layer-1)])

						s2 = sum(siteT3[1,I])

						ax4.plot([layer-1,layer],[s1,s2]*sf, color=col,zorder=5,alpha=0.7)

					end



					plot_arrows(ax1, dev_at, siteT; color=col, zorder=20,inds=I,min_abs_len=min_abs_len)
	
					plot_arrows(ax3, dev_at, siteT3; color=col, zorder=20,inds=I,min_abs_len=min_abs_len)
	
				end 
	
				ax4.set_xlim(0,nr_layers+1)

			end 	
					


#				@show sum(siteT[1,i_mid_slice]) sum(siteT3[1,i_mid_slice])

			println() 

			end 
	
			[a.set_ylim(0,1.1) for a in [ax2,ax4]]


			pad = 0.0

			[a.set_xlim(mu_min - pad*(mu_max-mu_min), mu_max + pad*(mu_max-mu_min)) for a in [ax2]]




			sleep(0.001)
			
			print("\r",Int(round(imu/length(MUS)*100)),"%  mu=$mu         ")

		end 
		println()	
	
		[a.legend() for a in [ax2]]
	

#Ax2_[2,2].legend()

end  

	
PyPlot.close(101) 
PyPlot.close(102) 

figsize = [21,12]*0.65 


fig1,Ax1 = PyPlot.subplots(2,4,num=101,figsize=figsize)#,sharex="col",sharey="col")

fig1.tight_layout()

sleep(0.1)

figure101(Ax1[1,:],Ax1[2,:])

