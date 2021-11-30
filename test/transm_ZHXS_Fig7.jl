include("transm_ZHXS.jl")
PyPlot.close(7)

fig,Ax = PyPlot.subplots(2,2,num=7,figsize=(10,7),sharex=true,sharey=true)

sleep(0.01)



function figure7(Ax2_)

	nr_at = 30#100
	nr_layers = 20 

	LAR = slices_ribbon(nr_layers, nr_at) 
	
	lead_latts,lead_labels,lead_conts = Utils.zipmap(two_lead_args(;LAR...)) do args 
	
		ll = square_latt_lead(args...; LAR...)
	
	
		return ll,args[3], lead_contacts(args[1], ll;  LAR...)
	
	end  .|> collect
	
	
	
	m = -0.5 
	delta = 1e-4


	BasisInfo = HoppingTerms.basis_info(HamiltBasis(true,true))
	lead_hopping = ZHXS_hopping(0, m, 0) 
	
	
	leads_ret = prep_lead.(lead_labels, lead_latts, [lead_hopping], delta) 
	leads_adv = prep_lead.(lead_labels, lead_latts, [lead_hopping], -delta)  
	



	el_proj, hole_proj = [Operators.Projector(BasisInfo(x), :atoms; nr_orb=length(BasisInfo(x)), nr_at=nr_at, dim=2) for x in "EH"]
	

	A, B = lead_labels 

	Deltas = LinRange(-1,1,30) 

	MUS =[0.2, 1.2, 1.4, 1.6]

	for (i_muS,(muS,ax2)) in enumerate(zip(MUS, Ax2_))


		println() 


		title = string("(",'a'+i_muS-1,") muS=$muS") 

		ax2.set_title(title)

		println(title)

		T,Tl,Tc = [zeros(length(Deltas)) for i=1:3]
	
	
#		ax2.plot(fill(sqrt(m^2 -Delta^2),2),[0,1.1],color="green",ls="--",alpha=0.5)
#		ax2.plot(extrema(MUS), fill(1/4,2),color="red",ls="--",alpha=0.5)
	
	
		for (i_D,Delta) in enumerate(Deltas)
	
			dev_hopping = ZHXS_hopping(Delta, m, muS) 
		
			NG_ret = NewGeometry(LAR, lead_conts, leads_ret, dev_hopping) 

			G_ret = GF(dev_hopping, NG_ret)(0) 

			se = get_SE(G_ret, NG_ret)
		

				NG_adv = NewGeometry(LAR, lead_conts, leads_adv, dev_hopping)
			
				G_adv = GF(dev_hopping, NG_adv)(0)
		
#				@assert G_ret(A,1,B,1)' ≈ G_adv(B,1,A,1)
		
	
		
	
			T[i_D] = ObservablesFromGF.CaroliConductance(G_ret, se, A, B, el_proj, el_proj)
	
			Tl[i_D] = real(ObservablesFromGF.CaroliConductance(G_ret, G_adv, se, A, A, el_proj, hole_proj))
	
	
			Tc[i_D] = ObservablesFromGF.CaroliConductance(G_ret, se, A, B, el_proj, hole_proj)

			for (t,lab,col) in [(T,"T","blue"),
															(Tc,"T_C","k"),(Tl,"T_L","red")]

				length(Deltas)<=50&& ax2.scatter(Delta,t[i_D],color=col, alpha=0.8,s=4)

				if i_D>1 
	
					I = i_D-1:i_D  

					kw = Dict{Symbol,Any}(:color=>col,:alpha=>0.8)

					i_D==2 && setindex!(kw, lab, :label)

					ax2.plot(Deltas[I], t[I]; kw...)

				end 
			end 
	
	
			ax2.set_ylim(0,1.1) 
			ax2.set_xlim(extrema(Deltas))
	
			sleep(0.001)
			
			print("\r",Int(round(i_D/length(Deltas)*100)),"%           ")

		end 
		println()	
	
	end	
	
	Ax2_[1,1].legend()

	[ax.set_xlabel("Delta") for ax in Ax2_[2,:]]

end  



figure7(Ax)
