using myLibs: Lattices, LayeredSystem, Algebra,GreensFunctions, Graph,ArrayOps,TBmodel,Utils 

using LinearAlgebra, SparseArrays 
import PyPlot


function get_hopp(nr_orb_)

	intercell = rand(ComplexF64,nr_orb_,nr_orb_)
#(ones(nr_orb_,nr_orb_)-ArrayOps.UnitMatrix(nr_orb_))
	function hopp(ri,rj)

		local_pot = ArrayOps.UnitMatrix(nr_orb_)*isapprox(norm(ri-rj),0,atol=1e-8)*2.0  
	
		atom_hopp = intercell*isapprox(norm(ri-rj), 1, atol=1e-8) 
	
		return local_pot + atom_hopp 
	
	end  

end 

function get_latt(Nxy::Vector{Int})

	l1 = Lattices.SquareLattice()

	Lattices.Superlattice!(l1, Nxy, Labels=first) 

end 

function get_atoms(Nxy::Vector{Int})

	l2 = get_latt(Nxy)

	return Lattices.PosAtoms(l2)  

end 

get_leadlatt(n::Int=1, label...) = Lattices.KeepDim!(
																					Lattices.Superlattice( 
																					 Lattices.SquareLattice(label...),
																					 [1,n]
																					 )
																					 , 1)
get_lead(L::Lattices.Lattice, args...) = GreensFunctions.PrepareLead(only(Lattices.sublatt_labels(L)), L, args...)

get_lead(n::Int, label::String, args...) = GreensFunctions.PrepareLead(label,get_leadlatt(n, label), args...)

d_nn = norm(diff(get_atoms([2,1]),dims=2))


isBond = Algebra.EuclDistEquals(d_nn; dim=2)

function default_LayerAtomRels(DL::Lattices.Lattice)

	atoms = Lattices.PosAtoms(DL)
	layers = parse.(Int, Lattices.labelToComponent(DL)) .+ 1

		# layer of each atom [1, 2, ..., nr_at] => [L(1), L(2), ..., L(nr_at)]
	
	la,ial = Utils.FindPartners(enumerate(layers), sortfirst=true)

	return Dict(

			:NrLayers=> maximum(layers),

			:LayerOfAtom => la,

			:IndsAtomsOfLayer => ial,

			:AtomsOfLayer => function AtomsOfLayer(L::Int)::AbstractMatrix
			
									Lattices.Vecs(atoms, ial(L))

										end
							)
end



@testset "G(g(data)) vs. G(g,data)" begin   

	for nr_orb in 1:3 
		
		for Nx=10:15:30, Ny=7:17:30 
		
			hopp = Dict(:Hopping=>get_hopp(nr_orb), :nr_orb=>nr_orb) 
	
			lead_width = max(div(min(Nx,Ny),2),1)
	
			atoms = get_atoms([Nx,Ny]) 
			
			all_bonds = mapreduce(hcat,findall(LinearAlgebra.triu!(isBond(atoms,atoms),1))) do C
	       [C.I[1],C.I[2]]
	      end
				
			for LAR0 in ["forced", default_LayerAtomRels(get_latt([Nx,Ny]))]
			
	
	
	
	#	leadA = get_lead(Lattices.Align_toAtoms!(get_leadlatt(lead_width,"A"),atoms,-1),hopp,0.001)
	
	#	leadB = get_lead(Lattices.Align_toAtoms!(get_leadlatt(lead_width,"B"),atoms,1),hopp,0.001)
		
				D,Slicer,LeadRels, = LayeredSystem.NewGeometry(
								atoms, LAR0; isBond=isBond, dim=2, nr_orb=nr_orb)
			
				nr_layers = D[:NrLayers]
			
				println() 
		
				@show (nr_orb,Nx,Ny,nr_layers,LAR0 isa String)
		
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
			
				g_noH = LayeredSystem.LayeredSystem_toGraph(nr_layers) # no H 
			
				g_withH = LayeredSystem.LayeredSystem_toGraph(HoppMatr, nr_layers)
			
			
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
		
					@test B1_2==B2_2 
		
					B2 = hcat(B2_1,B2_2)
		
		
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
			
				E1 = rand() -0.5 + 0.02im 
			
				G_old = GreensFunctions.GF_Decimation_fromGraph(E1, g_withH, Slicer) 
			
				G_new = GreensFunctions.GF_Decimation_fromGraph(E1, 
																												g_noH, 
																												data_H,
																												Slicer,
																												)
			
				GF_call_args = vcat( 
						[("Atom",i) for i in axes(atoms,2)],
						[("Layer",i) for i in 1:nr_layers]
					 )
			
				for a=GF_call_args, b=GF_call_args  
				
					@test G_old(a,b)≈G_old(a...,b...)≈G_old((a,b))
					@test G_old(a,b)≈G_new(a,b)≈G_new(a...,b...)≈G_new((a,b))
			
				end 
		
			end  
break 
		end 
#break 
	end 

end 

#LayeredSystem.Plot_Graph(("test10",D[:IndsAtomsOfLayer]) ,D[:NrLayers],g)



println()
error()

println()

@testset "lead basics" begin 
	Lead = get_leadlatt(1,"Lead-A")
	
	label = "A" 
	
	@test size(Lattices.LattVec(Lead))==(2,1)
	
	@test size.(GreensFunctions.PrepareLead(label, Lead)[:head])==[(2,1)]
	
	
	
	bridge = rand(2,3)
	
	@test size.(GreensFunctions.PrepareLead(label, Lead, bridge)[:head])==[(2,3),(2,1)]


	hopp(x,y=x) = rand(2,2)
	
	@test size.(GreensFunctions.PrepareLead(label, Lead, hopp, hopp, hopp)[:intracell])==[(2,2)]

	@test size.(GreensFunctions.PrepareLead(label, Lead, bridge, hopp, hopp, hopp)[:GF](rand())) == [(2,2),(2,2)]


println() 
end 



@testset "layered system geometry" begin 

	Nxy = [10,10]

	lead_width = max(div(minimum(Nxy),2),1)

	atoms = get_atoms(Nxy)

leadA = get_lead(Lattices.Align_toAtoms!(get_leadlatt(lead_width,"A"),atoms,-1))
leadB = get_lead(Lattices.Align_toAtoms!(get_leadlatt(lead_width,"B"),atoms,1))

#	PyPlot.close()
#	PyPlot.scatter(Lattices.eachcomp(atoms)...,label="Device")
#	PyPlot.scatter(Lattices.eachcomp(only(leadA[:head]))...,label="A")
#	PyPlot.scatter(Lattices.eachcomp(only(leadB[:head]))...,label="B")
#
#	PyPlot.legend()




	LayerAtomRels,LeadRels,VirtLeads = LayeredSystem.NewGeometry(atoms, "forced";
																														 isBond=isBond,
																														 Leads=[leadA,leadB],
																														 dim=2)

	nr_layers = LayerAtomRels[:NrLayers]



	g1 = LayeredSystem.LayeredSystem_toGraph(nr_layers, VirtLeads)

	@show nr_layers g1 

end 


@testset "layered system hopping" begin 

	Nxy = [9,5]
	nr_orb = 2  




	atoms = get_atoms(Nxy)
	lead_width = max(div(minimum(Nxy),2),1)


	hopp = Dict(:Hopping=>get_hopp(nr_orb), :nr_orb=>nr_orb)

	leadA = get_lead(Lattices.Align_toAtoms!(get_leadlatt(lead_width,"A"),atoms,-1),hopp,0.001)

	leadB = get_lead(Lattices.Align_toAtoms!(get_leadlatt(lead_width,"B"),atoms,1),hopp,0.001)

	@show keys(leadA) keys(leadB)



	LayerAtomRels,Slicer,LeadRels,VirtLeads = LayeredSystem.NewGeometry(
					atoms, "forced"; isBond=isBond, Leads=[leadA,leadB], dim=2,
					nr_orb=nr_orb)


	nr_layers = LayerAtomRels[:NrLayers] 

	HoppMatr(args...) = TBmodel.HoppingMatrix(LayerAtomRels[:AtomsOfLayer].(args)...; hopp...)

	g2 = LayeredSystem.LayeredSystem_toGraph(HoppMatr, nr_layers, VirtLeads)

	for layer=1:nr_layers

		@show size(LayeredSystem.get_graphH(g2, "Layer", layer))
		
		layer==1 && continue 
		
		@show size(LayeredSystem.get_graphH(g2, "Layer", layer-1, "Layer", layer))


	end 








end 




































nothing 
