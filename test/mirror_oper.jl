import myLibs:Operators ,Utils,TBmodel, ArrayOps, Lattices
import LinearAlgebra, PyPlot
using BenchmarkTools
using myLibs.Lattices: eachvec, eachcomp


nr_at_x = rand(4:10)

nr_at_y = rand(setdiff(4:10,[nr_at_x]))

nr_orb = rand(1:3)


mirror_plane = rand(1:2) 


nr_at = nr_at_x*nr_at_y 

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



#sort(rand(nr_wf))*2 .- 1 



dim = Lattices.VECTOR_STORE_DIM 
dim2 = Lattices.VECTOR_AXIS 


P = ArrayOps.init_zeros(ComplexF64, dim=>nr_wf, dim2=>size_H)

for i in 1:nr_wf

	selectdim(P,dim,i) .= LinearAlgebra.normalize(rand(ComplexF64,size_H))

end 


@testset "WF" begin 

	for psi in eachvec(P)

		@test LinearAlgebra.norm(psi)≈1 

	end 

end 

println() 

kwargs = (dim=dim,nr_at=nr_at,nr_orb=nr_orb,size_H=size_H)

@testset "op" begin 


	Op = Operators.Operator(M; kwargs...)


	@test Op.data ≈ kron(M,ArrayOps.UnitMatrix(nr_orb))


	@test Op.data^2 ≈ ArrayOps.UnitMatrix(size_H)

	@test  eltype(Op(P))<:Real


	Op2 = Operators.Mirror(mirror_plane, atoms; dim=dim, kwargs...)

	@test Op.data ≈ Op2.data 

	@test Op(P)≈ Op2(P) 










































end

