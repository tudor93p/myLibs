import myLibs: TBmodel, Operators, Lattices, H_Superconductor, ArrayOps
using LinearAlgebra, SparseArrays, BenchmarkTools 






function get_hopp(nr_orb_)

	function hopp(ri,rj)

		local_pot = ArrayOps.UnitMatrix(nr_orb_)*isapprox(norm(ri-rj),0,atol=1e-8)*2.0  
	
		atom_hopp = (ones(nr_orb_,nr_orb_)-ArrayOps.UnitMatrix(nr_orb_))*isapprox(norm(ri-rj), 1, atol=1e-8)*1.0
	
		return local_pot + atom_hopp 
	
	end  

end 


function get_PosAtoms(dim_,nr_at_,nr_orb_) 

	atoms = rand(dim_,max(1,div(nr_at_,dim_+1))); 

	return hcat(atoms, (atoms .+ setindex!(zeros(dim_),1,i) for i=1:dim_)...) 
	# such that neighbors exist 

end 

@testset "HoppingMatrix old==new" begin 

	times1::Float64 = 0 
	times2::Float64 = 0
	times3::Float64 = 0

	mem1::Float64 = 0 
	mem2::Float64 = 0 
	mem3::Float64 = 0 

	nr_at_min = 1000
	nr_at_max = Int(ceil(nr_at_min*1.)) 
	nr_orb_range = [2] 
	dim_range = [2] 

	for nr_at1=nr_at_min:nr_at_max 
		@show nr_at1 

		for dim=dim_range,nr_orb=nr_orb_range,nr_at2=nr_at1:nr_at_max
		
			#@show (dim,nr_orb,nr_at1,nr_at2)
		
			atoms_j = get_PosAtoms(dim,nr_at2,nr_orb) 
			atoms_i = get_PosAtoms(dim,nr_at1,nr_orb)
		
			kw = (Hopping=get_hopp(nr_orb), nr_orb=nr_orb)
			
			@assert size(atoms_i)!=size(atoms_j) || !isapprox(atoms_i,atoms_j,atol=1e-8)
		
			args = ((atoms_i,),(atoms_j,),(atoms_i,atoms_j),(atoms_j,atoms_i))
	
			H1 = [TBmodel.HoppingMatrix(a...; kw...) for a=args]
			H2 = [TBmodel.HoppingMatrix2(a...; kw...) for a=args]
			H3 = [TBmodel.HoppingMatrixAndNZ(a...; kw...) for a=args]
			
			times1 += @elapsed begin 
				[TBmodel.HoppingMatrix(a...; kw...) for a=args]
			end  
	
			times2 += @elapsed begin 
				[TBmodel.HoppingMatrix2(a...; kw...) for a=args]
			end  
			times3 += @elapsed begin 
				[TBmodel.HoppingMatrixAndNZ(a...; kw...) for a=args]
			end  

			mem1 += @allocated begin 
				[TBmodel.HoppingMatrix(a...; kw...) for a=args]
			end  
	
			mem2 += @allocated begin 
				[TBmodel.HoppingMatrix2(a...; kw...) for a=args]
			end 
			
			mem3 += @allocated begin 
				[TBmodel.HoppingMatrixAndNZ(a...; kw...) for a=args]
			end  

			for (h1,h2,(h3,nz3)) in zip(H1,H2,H3)
				
				@test h1≈h2 atol=1e-8 #skip=true
				@test h1≈h3 atol=1e-8 #skip=true 


				for (atom_i,atom_j,) in zip(SparseArrays.findnz(nz3)...)

					h3_ij = view(h3, 
											 TBmodel.Hamilt_indices(nr_orb,atom_i),
											 TBmodel.Hamilt_indices(nr_orb,atom_j),
											 )

					@test any(>(1e-6)∘abs,h3_ij)

				end  

				for (i,j,) in zip(SparseArrays.findnz(h3)...)

					@test nz3[TBmodel.orbat_indices(nr_orb, i)[2],
										TBmodel.orbat_indices(nr_orb, j)[2]
										]

				end 
 


	
			end 
	
		end    

		println()
		mem1_GB = mem1/=1e9
		mem2_GB = mem2/=1e9
		mem3_GB = mem3/=1e9
	
		@show times1 times2 times3 mem1_GB mem2_GB mem3_GB
	
		println()
	end 

end 



#
#TBL = ([0,1,2], [[1,0.],[2.,0],[0,1]], R)
#
#intra, (ms,Rms,Tms) = TBmodel.Compute_Hopping_Matrices(TBL,Hopping=hopp,
#																											 nr_orb=nr_orb)
#
