import myLibs: Lattices 
import Combinatorics

function gives_error(f::Function, args...; kwargs...)::Bool
	
	try 

		f(args...; kwargs...)

		return false 

	catch 

		return true 

	end 

end 

eq(x,y) = all(isapprox.(x,y,atol=1e-5))


@testset "Intact Lattice" begin 

	for space_dimension in 0:3, d in Combinatorics.powerset(1:space_dimension)
		

		alldims = 1:space_dimension

		lv = round.(rand(space_dimension,space_dimension),digits=3)
		
		at = round.(rand(space_dimension, 9),digits=3)
		
		L = Lattices.Lattice(lv, at)
		
		
		@test eq(L.LattDims, alldims)


		@test Lattices.isintact(L) && Lattices.dimsordered(L) 
		
		
		for mode in [:free, :absolute, :relative, :stored]

			@test Lattices.LattDim(L,mode)==space_dimension  
	
			@test eq(Lattices.LattDims(L,mode),alldims)
			
	
			@test gives_error(Lattices.LattDims, L, space_dimension+1, mode)
	
	
			@test eq(Lattices.LattDims(L, d, mode), d)
	
			@test eq(alldims,sort(union(d,Lattices.LattDims(L, d, mode; complement=true))))
	
	
			@test eq(d,	Lattices.LattDims(L, d, mode, :absolute))
	
	
			if !isempty(alldims) 
	
				randint = rand(alldims)
		
				@test eq(Lattices.LattDims(x->x[randint], L), L.LattDims[randint])
	
				if !isempty(d) 
	
					randint = rand(axes(d,1))
		
					@test eq(Lattices.LattDims(x->x[randint], L, d, mode) , d[randint])
	
				end 
	
			end 

		end 

	end 

end 
