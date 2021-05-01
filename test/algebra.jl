using myLibs: Algebra,Utils 



#A =rand(5)
#
#A[4]=A[2]
#
#foreach(println, Utils.Unique(A, inds=:all))
#foreach(println, Utils.Unique(A, dim=1,inds=:all))
#
#
#
#B = rand(3,5)
#
#B[:,4] = B[:,2]
#
#
#B[1,:] = B[2,:]
#
#println()
#
# 
#foreach(println, Utils.Unique(B, dim=2, inds=:all)[2])
#foreach(println, Utils.Unique(B, dim=1, inds=:all)[2])
#
#
#println()




#@testset "Lattices" begin 

#	latt = Lattices.SquareLattice()
#
#	@show size(latt.Sublattices["A"] )
#	@show Lattices.BondDirs(latt,1)

#end 

#for dim in [1,2]
for dim in []

	
	global dim2=[2,1][dim]

	global atoms = convert(Array{Float64},cat(cat.(collect(Base.product(1:3,1:5))[:]...,dims=dim)...,dims=dim2))
	
	@testset "get Bonds from atoms, dim=$dim" begin 
	
		for j=1:3
	
			b1 = Algebra.get_Bonds(atoms, 1, dim=dim, N=rand(axes(atoms,dim)))
			a1 = Algebra.bondRs_fromInds(b1, atoms, dim=dim)
	
			b2 = Algebra.get_Bonds(atoms, 1, dim=dim)
															

			a2 = Algebra.bondRs_fromInds(b2, atoms, dim=dim)


			b3 = Algebra.get_Bonds(atoms |> transpose, 1, dim=dim2, N=rand(axes(atoms,dim)))
	
			a3 = Algebra.bondRs_fromInds(b3, atoms|>transpose, dim=dim2)
			

			b4 = Algebra.get_Bonds(atoms |>transpose, 1, dim=dim2)

			a4 = Algebra.bondRs_fromInds(b4, atoms|>transpose, dim=dim2)

			@test b1==b2==b3==b4
	
			A2 = cat(cat.(a2...,dims=dim2)...,dims=dim) 
			A1 = cat(cat.(a1...,dims=dim2)...,dims=dim) 
			A3 = cat(cat.(a3...,dims=dim2)...,dims=dim) 
			A4 = cat(cat.(a4...,dims=dim2)...,dims=dim) 
	
			@test isapprox(A1,A2) & isapprox(A1,A3) & isapprox(A1,A4)
	
	
		end 
	
	end


	@testset "Bonds to matrix dim=$dim" begin 

		for j=1:3

			b1 = Algebra.get_Bonds(atoms, 1, dim=dim)

			a1 = Algebra.bondRs_fromInds(b1, atoms,dim=dim)

			b2 = Algebra.get_Bonds(atoms |>transpose, 1, dim=dim2)
			
			a2 = Algebra.bondRs_fromInds(b2, atoms|>transpose, dim=dim2)
			


			mb1 = Algebra.get_Bonds_toMatrix(b1; inds=true, dim=dim) 
			ma1 = Algebra.get_Bonds_toMatrix(a1; pos=true, dim=dim) 

			mb2 = Algebra.get_Bonds_toMatrix(b2; inds=true, dim=dim2)
			ma2 = Algebra.get_Bonds_toMatrix(a2; pos=true, dim=dim2) 

			mb1_ = Algebra.get_Bonds_toMatrix(b1, dim=dim, inds=true)
		
			ma1_ = Algebra.get_Bonds_toMatrix(atoms, mb1_, dim=dim, pos=true) 



			@test isapprox(mb1,transpose(mb2)) & isapprox(ma1,transpose(ma2)) &  isapprox(mb1,mb1_)

		end 


	end 


	@testset "matrix to Bonds dim=$dim" begin 

		for j=1:3

			b1 = Algebra.get_Bonds(atoms, 1, dim=dim)

			a1 = Algebra.bondRs_fromInds(b1, atoms;dim=dim)
																
#			b2,a2 = Algebra.get_Bonds(atoms |>transpose, 1, dim=dim2, pos=true)

			mb1 = Algebra.get_Bonds_toMatrix(b1; inds=true, dim=dim) 
			ma1 = Algebra.get_Bonds_toMatrix(a1; pos=true, dim=dim) 

			b1_ = Algebra.get_Bonds_fromMatrix(mb1, dim=dim,inds=true)

			a1_ = Algebra.get_Bonds_fromMatrix(ma1, dim=dim,pos=true)
			@test b1==b1_ && a1==a1_



		end 

	end 


	@testset "Bonds btw diff atoms dim=$dim" begin 

		for j=1:3

			b1  = Algebra.get_Bonds(atoms, 1, dim=dim)
			b2= Algebra.get_Bonds(atoms, atoms, 1, dim=dim)#,pos=true)
			b3= Algebra.get_Bonds(atoms, atoms, 1, dim=dim,order_pairs=true)

			@test b1==filter(x->x[2]>x[1],unique(b2))==b3


		end 

	end 
end





























































































































nothing
