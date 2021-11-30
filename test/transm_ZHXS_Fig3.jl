include("transm_ZHXS.jl")
PyPlot.close(3)
PyPlot.close(31)

fig2,Ax2 = PyPlot.subplots(2,2,num=3,figsize=(10,7),sharex=true,sharey=true)
fig1,Ax1 = PyPlot.subplots(2,2,num=31,figsize=(10,7),sharex=true,sharey=true)


sleep(0.1)



function figure3a(Ax1_,Ax2_)

	nr_layers = 20 

	widths = [80,100,150,100] 

	widths = div.(widths,5)

	BasisInfo = HoppingTerms.basis_info(HamiltBasis(true,true))

	MUS = LinRange(0,2,30)


	for (i_w,(nr_at,ax1, ax2)) in enumerate(zip(widths, Ax1_, Ax2_))


		println() 

		@show nr_at 



		LAR = slices_ribbon(nr_layers, nr_at) 
	
#		nr_at==widths[1] && 
		LayeredSystem.plot_layers(ax1; LAR..., s=3)
	
		lead_latts,lead_labels,lead_conts = Utils.zipmap(two_lead_args(;LAR...)) do args 
	
			ll = square_latt_lead(args...; LAR...)
	
			for n = 0:1
	
				Lattices.plot_atoms(Lattices.Atoms_ManyUCs(ll; Ns=n),ax1;s=3,max_atoms=40)
	
			end  
	
	
			return ll,args[3], lead_contacts(args[1], ll;  LAR...)
		
		end  .|> collect
	
	
	
		delta = 1e-4

		Delta = 0.35  

		m = i_w==4 ? -0.8 : -0.5 

	

		title = string("(",'a'+i_w-1,") $nr_layers x $nr_at")


		ax1.set_title(title)
		ax2.set_title("$title, m=$m")
			

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
		
#				@assert G_ret(A,1,B,1)' ≈ G_adv(B,1,A,1)
		
	
		
	
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
	
	
			ax2.set_ylim(0,1.1) 
			ax2.set_xlim(extrema(MUS))
	
			sleep(0.001)
			
			print("\r",Int(round(imu/length(MUS)*100)),"%           ")

		end 
		println()	
	
	end	
	
	Ax2_[2,2].legend()

end  



figure3a(Ax1, Ax2)
