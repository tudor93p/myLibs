using myLibs: Algebra,Lattices ,Utils 


@testset "Superlattice" begin 


	latt = Lattices.SquareLattice()


	@show size(Lattices.PosAtoms(latt))

	@show Lattices.sublatt_labels(latt)


#	n = rand(1:21,2,2)

n = [3,5]

	l1 = Lattices.Superlattice(latt,n)#[5,1])

	@show size(Lattices.PosAtoms(l1))
	@show Lattices.LattVec(l1)
	


	@show Lattices.sublatt_labels(l1)


	println()
	println()

	l2 = Lattices.to_myODict(["A"=>rand(2,1),rand(2,5),Dict(2=>rand(2,3)),(["A","B"],rand(2,2)),rand(2,4)],l1)

@show	keys(l2)
	function test(A,B)

		R,P = Lattices.PosAtoms(A), Lattices.PosAtoms(B)
		
		@test size(R)==size(P)

		if size(R)==size(P) 

		d= Algebra.OuterDist(R,P,dim=2)

		@test all(count.(eachcol(isapprox.(d,0))).==1)

	end 

		@test isapprox(A.LattVec,B.LattVec)
	end 

	@show size.(values(l2),2)


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


	Lattices.ShiftAtoms(l5, n=n)
	Lattices.ShiftAtoms!(l5, n=n)



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


Lattices.ReduceDim!(l5,1)

@show Lattices.ReciprocalVectors(l5)

@show Lattices.BrillouinZone(l5)


println()
println()





n = [3,5]

L = Lattices.Superlattice(Lattices.SquareLattice("A"), [3,2], Labels=x->x[1])


@show L 



@show Lattices.labelToComponent(L, labels_contain="B")










































































































































nothing
