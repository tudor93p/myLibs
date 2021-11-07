import myLibs: Lattices 
import Combinatorics, Random 

function gives_error(f::Function, args...; kwargs...)::Bool
	
	try 

		f(args...; kwargs...)

		return false 

	catch 

		return true 

	end 

end 

eq(x::Union{AbstractArray,Number},y::Union{AbstractArray,Number})::Bool = all(isapprox.(x,y,atol=1e-5))

function eq(L1::Lattices.Lattice, L2::Lattices.Lattice)::Bool

	eq(Lattices.PosAtoms(L1), Lattices.PosAtoms(L2)) && eq(L1.LattVec, L2.LattVec) && eq(L1.LattDims, L2.LattDims) 

end 

@testset "Intact Lattice" begin 

	for space_dimension in 0:3 
		
		alldims = 1:space_dimension

		lv = round.(rand(space_dimension,space_dimension),digits=3)
		
		at = round.(rand(space_dimension, 9),digits=3)
		
		L = Lattices.Lattice(lv, at)
		
		
		@test eq(L.LattDims, alldims)

		@test Lattices.isintact(L) && Lattices.dimsordered(L) 
	

		L0 = Lattices.Lattice(lv, at) 

		L1 = Lattices.ReduceDim(Lattices.Lattice(lv, at))


		@test Lattices.LattDim(L1)==0 

		@test space_dimension==Lattices.NrVecs(L0.LattVec)

		Lattices.ReduceDim!(L0) 

		@test eq(L0,L1)


		
		for mode in [:free, :absolute, :relative, :stored],d in Combinatorics.powerset(alldims)

			@test Lattices.LattDim(L,mode)==space_dimension  
	
			@test eq(Lattices.LattDims(L,mode),alldims)
			
	
			@test gives_error(Lattices.LattDims, L, space_dimension+1, mode)
	
	
			@test eq(Lattices.LattDims(L, d, mode), d)
	
			@test eq(alldims,sort(union(d,Lattices.LattDims(L, d, mode; complement=true))))
	
	
			@test eq(d,	Lattices.LattDims(L, d, mode, :absolute))
	
	
			isempty(alldims) && continue 
	
			randint = rand(alldims)
		
			@test eq(Lattices.LattDims(x->x[randint], L), L.LattDims[randint])
	
			isempty(d) && continue 
	
			randint = rand(axes(d,1))
		
			@test eq(Lattices.LattDims(x->x[randint], L, d, mode) , d[randint])



		end 

	end 

end  



@testset "Reduced Lattice" begin 

	for space_dimension in [2,3]

		alldims = 1:space_dimension
		
		lv = round.(rand(space_dimension,space_dimension),digits=3)
		
		at = round.(rand(space_dimension, 9),digits=3)

		L1 = Lattices.Lattice(lv, at)

		L2 = Lattices.ReduceDim(L1, 2, :relative)


		Lattices.ReduceDim!(L1, 2, :relative)

		@test eq(L1.LattDims, vcat(1,3:space_dimension))

		@test eq(L1,L2)



		L3 = Lattices.ReduceDim(L1, 1, :relative) 
		
		Lattices.ReduceDim!(L1, 1, :relative) 

		@test eq(L1,L3) 

		@test eq(L1.LattDims, 3:space_dimension)
		

		@test gives_error(Lattices.ReduceDim, L1, 1, :absolute)


		for L in [L1,L2,L3], i in Combinatorics.powerset(alldims)

			for j in Random.shuffle(collect(Combinatorics.powerset(alldims))), m in Random.shuffle([:relative, :absolute, :stored])


				Lr = Lattices.ReduceDim(Lattices.Lattice(lv, at), j, m)

				@test eq(Lattices.ReduceDim(L, i, :stored), #erases history
								 Lattices.ReduceDim(Lr, i, :stored)
								 )
				@test eq(Lattices.ReduceDim!(L, i, :stored), #erases history
								 Lattices.ReduceDim(Lr, i, :stored)
								 )
			end  


			@test eq(Lattices.Lattice(lv,at),
							 Lattices.ReduceDim!(L, Int[], :stored)) # reset latts 


		end

		Lattices.ReduceDim!(L1, 2, :relative) 
		Lattices.ReduceDim!(L2, 2, :relative)

		#x (and z) are left 


		for i in alldims 
			
			for m in [:relative, :absolute]

				test = (m==:relative && i<space_dimension) | (m==:absolute && i!=2)

				@test xor(test, gives_error(Lattices.KeepDim, L1, i, m))

				if test 

				end 
	
		end 

			@test !gives_error(Lattices.KeepDim, L1, i, :stored)


			Lattices.KeepDim!(L3, i, :stored)
			@test eq(L3,Lattices.KeepDim(L1, i, :stored))


		end 

		@test eq(Lattices.KeepDim(L1, 1, :relative), 
						 Lattices.KeepDim(L2, 1, :absolute))


		@test eq(Lattices.KeepDim(L1, 1:(space_dimension-1), :relative), 
						 Lattices.KeepDim(L2, vcat(1,3:space_dimension), :absolute))


		space_dimension>2 && @test eq(Lattices.KeepDim(L1, 2, :relative),
																	Lattices.KeepDim(L2, 3, :absolute))




#Lattices.KeepDim!()
		

	

	end 

end 
