import myLibs:Operators ,Utils,TBmodel, ArrayOps, Lattices,Algebra,HoppingTerms
import LinearAlgebra, PyPlot, Random ,SparseArrays
using BenchmarkTools
using myLibs.Lattices: eachvec, eachcomp
LA = LinearAlgebra

Random.seed!(Int(round(time()*100)))

nr_at_x = rand(1:10)

nr_at_y = rand(1:10) 

nr_at = nr_at_x*nr_at_y 

#nr_orb = rand(1:10)

nr_orb = rand([rand(1:10),nr_at + 1])

nr_orb = 2

@show nr_at nr_orb 

#@show nr_orb <= nr_at 

mirror_plane = rand(1:2) 


dim = Lattices.VECTOR_STORE_DIM 
dim2 = Lattices.VECTOR_AXIS 


size_H = nr_at*nr_orb 

nr_wf = Int(round(size_H*(0.2+rand()*0.6))) 

static_dirs = setdiff(1:2, [mirror_plane])


atoms = Lattices.VecsOnDim(rand(1,2) .+ hcat(
						Algebra.FlatOuterSum(1:nr_at_x,fill(0,nr_at_y)),
						Algebra.FlatOuterSum(fill(0,nr_at_x),1:nr_at_y));dim=1)


@assert Lattices.NrVecs(atoms)==nr_at

center  = x_mid,y_mid = Algebra.Mean(atoms,2)  

#atoms .-= center 



#PyPlot.close()  


#fig, ax1 = PyPlot.subplots(1)

#Lattices.plot_atoms(atoms, ax1) 
#ax1.plot(fill(x_mid,2),extrema(atoms_y),c="k",lw=1)
#atoms_x,atoms_y = eachcomp(atoms)  
#ax1.plot(extrema(atoms_x), fill(y_mid,2), c="k",lw=1)

mirror_position = center[mirror_plane] 

#mirrored_atoms = copy(atoms) 
#
#selectdim(mirrored_atoms, dim2, mirror_plane) .= (2mirror_positon .- selectdim(atoms, dim2, mirror_plane) )
#
#ax1.scatter(eachcomp(mirrored_atoms)...,s=1)
#

M = Lattices.MirrorReflectionMatrix(atoms, mirror_plane, mirror_position)

kwargs = (dim=dim,nr_at=nr_at,nr_orb=nr_orb,size_H=size_H)


false && @testset "Combine operators" begin

	for i=0:3 
	
#		@show i mirror_plane atoms kwargs Algebra.PauliMatrix(i) 

		Op0987 = Operators.SpecialMirror(i, mirror_plane, atoms; kwargs...)

		@test Op0987.data ≈ kron(M,Algebra.PauliMatrix(i))

	end 


	for M_ in [ArrayOps.UnitMatrix(nr_at),
#								 rand(Bool, nr_at, nr_at),
#								 rand(Float64, nr_at, nr_at),
								 rand(ComplexF64, nr_at, nr_at),
								 SparseArrays.sprand(ComplexF64, nr_at, nr_at,0.3),
								 rand(nr_at),
								 ], orb_op in [
															 fill(1,nr_orb),
															 ArrayOps.UnitMatrix(nr_orb),
															 SparseArrays.spdiagm(fill(1,nr_orb)),
#								 rand(Float64, nr_orb, nr_orb),
#								 rand(ComplexF64, nr_orb, nr_orb),
								 ]

 		a = Operators.repeatOperator_manyOrbitals(M_, nr_at, nr_orb, size_H)
#		@time Operators.repeatOperator_manyOrbitals(M_, nr_at, nr_orb, size_H)

		b = Operators.combineOperators_AtOrb(M_, orb_op)
#		@time Operators.combineOperators_AtOrb(M_, orb_op)


		if ndims(a)==ndims(b) 
			
			@test a≈b

		elseif ndims(a)==2 && ndims(b)==1 

			@test LA.diag(a)≈b 

		elseif ndims(a)==1 && ndims(b)==2

			@test a≈LA.diag(b)

		end 


		Op = Operators.SpecialMirror(Array(orb_op), mirror_plane, atoms,
																 mirror_position; kwargs...)



		d_ = Operators.combineOperators_AtOrb(M, orb_op)

		@test Op.data ≈ d_

		if Op.data ≈ d_ 
			
		else 

			@show size(d_) size(Op.data)

			@show M orb_op 

			break 

		end 


#		@show M_ isa SparseArrays.AbstractSparseMatrix 
#		@show orb_op isa SparseArrays.AbstractSparseMatrix 
#		@show b isa SparseArrays.AbstractSparseMatrix 

#		println() 

	end 



end 

#sort(rand(nr_wf))*2 .- 1 

println() 



P = ArrayOps.init_zeros(ComplexF64, dim=>nr_wf, dim2=>size_H)

for i in 1:nr_wf

	selectdim(P,dim,i) .= LinearAlgebra.normalize(rand(ComplexF64,size_H))

end 


false && @testset "WF" begin 

	for psi in eachvec(P)

		@test LinearAlgebra.norm(psi)≈1 

	end 

end 

println() 


false && @testset "op" begin 


	Op = Operators.Operator(M; kwargs...)


	@test Op.data ≈ kron(M,ArrayOps.UnitMatrix(nr_orb))


	@test Op.data^2 ≈ ArrayOps.UnitMatrix(size_H)

	@test  eltype(Op(P))<:Real


	Op2 = Operators.Mirror(mirror_plane, atoms; kwargs...)

	@test Op.data ≈ Op2.data 

	@test Op(P)≈ Op2(P) 


	for i=0:3 

		Op3 = Operators.SpecialMirror(i, mirror_plane, atoms;  kwargs...)

		@show i unique(LA.eigvals(Op3.data))

		@test isreal(Op3(P))

	end 









end 




@testset "Hamiltonian basis" begin 


	for use_spin in [false, true], use_Nambu in [false, true] 
		
		println("\n ====== spin: $use_spin  ======= Nambu: $use_Nambu ===== \n")





		hb = HoppingTerms.HamiltBasis(use_spin,use_Nambu)


		kwargs = (dim=dim,nr_at=nr_at,nr_orb = hb.matrix_dim,
							size_H = nr_at*hb.matrix_dim)

		@show hb.spin hb.charge  


		bi = HoppingTerms.basis_info(hb)

		@show bi("QP")
		@show bi("E")
		@show bi("H")

		if use_spin & use_Nambu 

#			kron(Nambu, spin)
			@test bi("E")≈kron([1,0],[1,1])

			@test hb.spin ≈ kron([1,1],[1,-1])


			@test repeat([1,-1], outer=2)≈kron([1,1],[1,-1])
			@test repeat([1,-1],inner=2)≈kron([1,-1],[1,1]) 

			println("Kronecker product atoms x Nambu x spin") 

		end 

		println()

		get_info = HoppingTerms.basis_info(hb)
#
#		@show HoppingTerms.representation(hb)
#		@show get_info()
#
#		println()
#
#		@show HoppingTerms.representation(hb; spin=rand(2))
#		@show HoppingTerms.representation(hb; spin=rand(2,2))
#		@show get_info(;spin=2)
#
#		println() 
#		
#		@show HoppingTerms.representation(hb; Nambu=rand(2))
#		@show HoppingTerms.representation(hb; Nambu=rand(2,2))
#		@show get_info(; Nambu=1)
#
#		println()
#
#		@show HoppingTerms.representation(hb; spin=rand(2,2), Nambu=rand(2)) 
#		@show HoppingTerms.representation(hb; spin=rand(2,2), Nambu=rand(2,2)) 
#		@show get_info(; Nambu=1, spin=2)  
#
#		println()


#@show kwargs 
		for i = 0:3 
			
#			@show i 

			PH = get_info(;Nambu=i)

#			@show size(PH )

			Op637 = Operators.Operator(PH, :atoms; kwargs...)

			@show Op637 == Operators.Operator(PH; kwargs...)

#			@show Op637.diag size(Op637.data) Op637.sums Op637.inds Op637.isHermitian maximum(abs, Op637.data) 
			
		end 



		println()
	end 











end 

