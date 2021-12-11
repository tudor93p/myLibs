using myLibs: Algebra, Lattices, Utils 

import LinearAlgebra
@testset "Superlattice" begin 


	latt = Lattices.SquareLattice()


	@show size(Lattices.PosAtoms(latt))

	@show Lattices.sublatt_labels(latt)

	@show Lattices.get_Bonds(latt)

#	n = rand(1:21,2,2)

n = [3,5]

	l1 = Lattices.Superlattice(latt, n)#[5,1])



	
	@show size(Lattices.PosAtoms(l1)) 

	@show Lattices.LattVec(l1)
	


	@show Lattices.sublatt_labels(l1)


	println()
	println()

	l2 = Lattices.to_myODict(["A"=>rand(2,1),
														rand(2,5),
														Dict(2=>rand(2,3)),
														(["A","B"],rand(2,2)),
														rand(2,4)
														],l1)

@show	keys(l2)
	function test(A,B)

		R,P = Lattices.PosAtoms(A), Lattices.PosAtoms(B)
		
		@test size(R)==size(P)

		if size(R)==size(P) 

			d = Lattices.OuterDist(R, P)

			@test all(count.(eachcol(isapprox.(d,0))).==1)

		end 

		@test isapprox(A.LattVec,B.LattVec)
	end 

	@show Lattices.NrVecs.(values(l2))


	println()
	println()

	get_label(n::AbstractVector{Int}) = n[2]#sum(n)#[2]

	@show n

	l3 = Lattices.Superlattice(latt,n, Labels=get_label)
	

	@show Lattices.sublatt_labels(l3)

	@show size.(values(l3.Atoms),2)



	test(l1,l3)

	println()
	println()



	l4 = Lattices.Superlattice([latt],[n])

	 test(l1,l4)

	l5 = Lattices.SquareLattice()

	Lattices.Superlattice!(l5, n, Labels=get_label)

	 test(l1, l5)

end 

@testset "ShiftAtoms" begin 

	l5 = Lattices.SquareLattice()

	n = [rand(1:10),rand(1:10)]

a = 	 Lattices.PosAtoms(l5)

b = Lattices.ShiftAtoms(l5, n=n) |> Lattices.PosAtoms

c = 	Lattices.ShiftAtoms!(l5, n=n) |> Lattices.PosAtoms

	@test isapprox(b,c)

#	Lattices.map(l5) do kind::Symbol 
#
#		length(string(kind))
#
#	end  |> println 
#
#	println()
#	
#
#
#
#	Lattices.map(l5) do D::AbstractMatrix 
#
#		length(D)
#
#	end  |> println 



end 

l5 = Lattices.SquareLattice()



for item in Lattices.NearbyUCs(l5)

	println()

	foreach(println, eachcol(item))

	println()

end

println("Reduce Dim")

@show Lattices.LattDims(l5) size(Lattices.LattVec(l5)) 

println.(eachcol(Lattices.LattVec(l5)))


Lattices.ReduceDim!(l5, 1)

@show Lattices.LattDims(l5,:absolute) 
@show Lattices.LattDims(l5,:relative)
@show Lattices.LattDims(l5,:stored)

@show size(Lattices.LattVec(l5))
println.(eachcol(Lattices.LattVec(l5)))

println()


@show Lattices.ReciprocalVectors(l5)


@testset "ReciprocalVectors" begin 

for (x,y) in zip(Lattices.eachvec(Lattices.LattVec(l5)),Lattices.eachvec(Lattices.ReciprocalVectors(l5)))

	@test isapprox(LinearAlgebra.dot(x,y),2pi)

end 
println()

@show Lattices.BrillouinZone(l5)


println()
end 

println()






n = [3,5]

L = Lattices.Superlattice(Lattices.SquareLattice("A"), [3,2], Labels=x->x[1])

for (k,v) in  L.Atoms 

	@show k v 

	println()

end 

println()



@show Lattices.labelToComponent(L, labels_contain="1") 
@show Lattices.PosAtoms(L, labels_contain="1") 



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

println()

L = Lattices.Superlattice(Lattices.SquareLattice("A"), [[4 3];[1 2]]) 

@show Lattices.LattVec(L)

println()

@show Lattices.LattDim(L)

@show Lattices.ReciprocalVectors(L) |> size 

for (a,b) in zip(eachcol(Lattices.ReciprocalVectors(L)), eachcol(Lattices.LattVec(L))) 

#	@show a b 
@test isapprox(2pi,LinearAlgebra.dot(a,b))
#			 println()
		 end 


println()
Lattices.KeepDim!(L, rand(1:2))

@show Lattices.LattDim(L) 


@show Lattices.ReciprocalVectors(L) |> size

@testset "ReciprocalVectors" begin 


for (a,b) in zip(eachcol(Lattices.ReciprocalVectors(L)), eachcol(Lattices.LattVec(L)))
#	@show a b 
@test isapprox(2pi,LinearAlgebra.dot(a,b))
#			 println()
       end


		 end 



































































































































nothing
