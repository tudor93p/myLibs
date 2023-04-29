using myLibs: Lattices, LayeredSystem, Algebra,GreensFunctions, Graph,ArrayOps,TBmodel,Utils 

using LinearAlgebra, SparseArrays 
#import PyPlot

include("mock_DevLeads.jl")

@testset "G(g(data)) vs. G(g,data)" begin   

	for nr_orb in [1,3]
		
		hopp = Dict(:Hopping=>get_hopp(nr_orb), :nr_orb=>nr_orb)  

		for Nx=10:15:20, Ny=7:17:20 
	
			atoms = get_atoms([Nx,Ny])  

			leads_  = get_two_leads(max(div(Ny,2),1), atoms; hopp...)

			
			all_bonds = mapreduce(hcat,findall(LinearAlgebra.triu!(isBond(atoms,atoms),1))) do C
	       [C.I[1],C.I[2]]
	      end
				
				for LAR0 in ["forced", default_LayerAtomRels(get_latt([Nx,Ny]))],Leads in [[],leads_]
			
				LAR0=="forced"&&!isempty(Leads)&& continue 	


				for L in Leads  
					for h in L[:intracell]
						@test ishermitian(h)
					end 
					for h in L[:intercell]
						@test !ishermitian(h)
					end 
				end  

				D,Slicer,LeadRels, VirtLeads = LayeredSystem.NewGeometry(
								atoms, LAR0; isBond=isBond, dim=2, nr_orb=nr_orb,
								Leads=Leads) 

				for (k,L) in VirtLeads  
					for h in L[:intracell]
						@test ishermitian(h)
					end 
					for h in L[:intercell]
						@test !ishermitian(h)
					end 
				end 
		

				nr_layers = D[:NrLayers]
			
				println() 
		
				@show (nr_orb,Nx,Ny,nr_layers,LAR0 isa String,length(VirtLeads))

		
				for i=axes(atoms,2)
			
					@test 1<=D[:LayerOfAtom](i)<=nr_layers
			
				end  
			
				for l=1:nr_layers 
			
					@test all(1 .<= D[:IndsAtomsOfLayer](l) .<= size(atoms,2))
			
					for j=1:nr_layers
						
						@test xor(l==j,isdisjoint(D[:IndsAtomsOfLayer](l),
																			D[:IndsAtomsOfLayer](j)
																			))
					end 
			
					D[:AtomsOfLayer](l)
				end 
			
			
			
			
				HoppMatr(args...) = TBmodel.HoppingMatrix(D[:AtomsOfLayer].(args)...; hopp...)

				g_noH = LayeredSystem.LayeredSystem_toGraph(nr_layers, VirtLeads) 

					for ll in GreensFunctions.get_leadlabels(g_noH)

						@test !isempty(VirtLeads)

						@test !isempty(LayeredSystem.get_graphH(g_noH, ll,1))
						

					end 

			
				g_withH = LayeredSystem.LayeredSystem_toGraph(HoppMatr, nr_layers, VirtLeads)
				
					for ll in GreensFunctions.get_leadlabels(g_withH)

						@test !isempty(VirtLeads)
						@test !isempty(LayeredSystem.get_graphH(g_withH, ll,1))

					end 

		
				leadlabels = Graph.get_prop(g_noH,:LeadLabels) 

				@test leadlabels==Graph.get_prop(g_withH,:LeadLabels)
			
				data_H,data_B = LayeredSystem.condStore_sharedHB(D, atoms; hopp...)
		
				@test !isempty(data_H) 
				@test !isempty(data_B)
		
				for l=2:nr_layers
					I = D[:IndsAtomsOfLayer](l-1)
					J = D[:IndsAtomsOfLayer](l)
		
					B1_1 = all_bonds[:,all(in(J),all_bonds,dims=1)[:]] 
					
					B1_2 = all_bonds[:,map(eachcol(all_bonds)) do (i,j)
												
													 in(i,I)&in(j,J)
												 end[:]]
					
		
					B1 = all_bonds[:,map(eachcol(all_bonds)) do (i,j)
											
		
													 in(j,J) & (in(i,I)|in(i,J))
		
												 end[:]]
				
					B2_1 = get(data_B, (l,l), zeros(Int,2,0)) 
					B2_2 = get(data_B, (l-1,l),zeros(Int,2,0))
		
					for (i2,b2)=enumerate(B2_1)
						B2_1[i2]=J[b2]
					end  
		
					@test B1_1==B2_1 
		
					for (j2,b2)=enumerate(eachcol(B2_2)) 
		
						B2_2[1,j2] = I[b2[1]]
						B2_2[2,j2] = J[b2[2]]
		
					end 
		
					@test size(B1_2)==size(B2_2)

					for b1 in eachcol(B1_2)
		
						@test any(==(b1), eachcol(B2_2))
		
					end

		
					B2 = hcat(B2_1,B2_2)

					@test ishermitian(get(data_H, (l,l), hcat(0)))
		
		#			@show size(get(data_H, (l,l), []))
		#			@show size(get(data_H, (l-1,l), []))
		
		#			println("B1:",eachcol(B1)...)
		#			println("B2:",eachcol(B2)...)
		
					@test size(B1)==size(B2)	 
		
					for b1 in eachcol(B1) 
		
						@test any(==(b1), eachcol(B2))
		
					end
		
				end 
		
		
				#		@show keys(data_B) keys(data_H)
				
			#-----  	 test GF ---- # 
			
				E1 = rand() -0.5 
			
#				E1 = 0.3232 

				if isempty(VirtLeads)
					E1 += 0.01im 
				end 
			
				G_old = GreensFunctions.GF_Decimation_fromGraph(E1, g_withH, Slicer;
																		leads_have_imag=true) 

				G_old2 = GreensFunctions.GF_Decimation(HoppMatr, nr_layers, VirtLeads;
																							 translate=Slicer,
																							 leads_have_imag=true,
																							 )(E1)

				G_old3 = GreensFunctions.GF_Decimation(hopp, VirtLeads, Slicer;
																							 D...,
																							 leads_have_imag=true
																							 )(E1)
			
				G_new = GreensFunctions.GF_Decimation_fromGraph(E1, 
																												g_noH, 
																												data_H,
																												Slicer;
																												leads_have_imag=true,
																												)
		
				G_new2 = GreensFunctions.GF_Decimation(data_H, nr_layers, VirtLeads;
																							 translate=Slicer,
																							 leads_have_imag=true,
																							 )(E1)

				GF_call_args = vcat( [(l,1) for l in leadlabels],
						[("Atom",i) for i in axes(atoms,2)],
						[("Layer",i) for i in 1:nr_layers]
					 ) 

			
				for a=GF_call_args, b=GF_call_args
				
						G0 = G_old(a,b)
					for ab=((a,b),(a...,b...)), fG=(G_old,G_old2,G_new,G_new2,G_old3)

						@test G0≈fG(ab...)
					
					end 
				end 
		
			end  
#break 
		end 
#break 
	end 

end 

#LayeredSystem.Plot_Graph(("test10",D[:IndsAtomsOfLayer]) ,D[:NrLayers],g)



println() 








































































































nothing 
