using myLibs: TBmodel, Operators, Lattices, H_Superconductor, ArrayOps
using LinearAlgebra 


R = rand(2,3);
R= hcat(R,R.+[0,1],R.+[1,0])#,R.+[0,-1],R.+[-1,0])

nr_orb = 2

local_pot(ri,rj) = ArrayOps.UnitMatrix(nr_orb)*isapprox(norm(ri-rj),0)*2.0  

atom_hopp(ri,rj) = (ones(nr_orb,nr_orb)-ArrayOps.UnitMatrix(nr_orb))*isapprox(norm(ri-rj), 1)*1.0

hopp(ri,rj) = local_pot(ri,rj) + atom_hopp(ri,rj) 


TBmodel.HoppingMatrix(R; Hopping=hopp, nr_orb=nr_orb)

TBL = ([0,1,2], [[1,0.],[2.,0],[0,1]], R)

intra, (ms,Rms,Tms) = TBmodel.Compute_Hopping_Matrices(TBL,Hopping=hopp,
																											 nr_orb=nr_orb)


@show size(intra)

@show ms 
@show Rms 

@show size.(Tms)



for (argH, k) in [("k", rand(2)), ("phi", rand())]

	println()

	@show argH 

	@show size(TBmodel.Bloch_Hamilt(TBL; Hopping=hopp, argH=argH,nr_orb=nr_orb)(k))


	dir = length(k)

	v = TBmodel.Assemble_Velocity_Operator((ms, Rms, Tms), dir; argH=argH)

	v2 = TBmodel.Bloch_Velocity(TBL, dir; Hopping=hopp, argH=argH,nr_orb=nr_orb)

	@show size(v(k)) 

	@assert isapprox(v(k),v2(k))

end 


println() 


H = TBmodel.Bloch_Hamilt(TBL; Hopping=hopp,nr_orb=nr_orb)



#@show H(1) - H([2]) + H([3,4]) - H([5,6,7])

@show H(1)

@show size(H(1))


println()

F = Operators.Position(1, R, nr_orb=nr_orb,dim=2) 


E =  LinearAlgebra.eigen(Matrix(H(1)))
												 

#@show propertynames(E)


P = E.vectors;


@assert isapprox(H(1)*P[:,1:1] ,E.values[1]*P[:,1:1])


@show size(F(P[:,1:3])) 


@show size(H(1))
#@show size(F(transpose(P[:,1:3]);dim=1))





#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#




println("\n\n------------ Parallel and perpedicular Hamilt ------\n")



#L = Lattices.KeepDim(Lattices.Superlattice(Lattices.SquareLattice(),[3,1]), 2)
L = Lattices.Superlattice(Lattices.SquareLattice(),[3,1]) 

#Lattices.KeepDim!(L, 2)



@show Lattices.PosAtoms(L) |>size
@show Lattices.LattDim(L) 
@show Lattices.VecLen(L)
@show Lattices.LattDims(L)

println() 

tbl = Lattices.NearbyUCs(L)

println.(tbl)


intra,inter = TBmodel.Compute_Hopping_Matrices(tbl; Hopping=hopp,nr_orb=nr_orb)

println()


@show intra 

println("\ninter:")
!isnothing(inter) && println.(inter)

println()



H = TBmodel.Bloch_Hamilt(tbl; Hopping=hopp,nr_orb=nr_orb)


@show H(rand(Lattices.LattDim(L)))


println()


Hpar = TBmodel.BlochHamilt_Parallel(L, Dict(:Hopping=>hopp,:nr_orb=>nr_orb))


@show Hpar(rand()) #Bloch phase #Lattices.LattDim(L)-1))

println()



Hperp = TBmodel.BlochHamilt_Perpendicular(L, Dict(:Hopping=>hopp,:nr_orb=>nr_orb))



@show Hperp[1]




HPar,Hperp = TBmodel.BlochHamilt_ParallPerp(L, Dict(:Hopping=>hopp,:nr_orb=>nr_orb))




@show Hperp[1]


@show Hpar(rand())#Bloch phase 



println() 
println() 

latt = Lattices.Lattice([1.0 0.0; 0.0 1.0], [0.0; 0.0], nothing, [2])

HParams = (Hopping = -1.0, ChemicalPotential = 2.0, SC_Gap = nothing) 

HParams = (Hopping = -1.0, ChemicalPotential = 2.0, SC_Gap = nothing,
					 LocalPotential= function v(r::AbstractVector{<:Real})
						 r[1]
						end 
					) 

hopp2 = H_Superconductor.SC_Domain(HParams, [1.0])

@testset "t(i,i+1)==1" begin 

	for r in eachcol(rand(2,5)) 

		@test hopp2[:Hopping](r,r)[1,1] ≈ r[1]+HParams.ChemicalPotential
	
		for i=1:2, s=[-1,1] 
	
			r1 = zeros(2) 
			r1[i] = s
	
			for args in [(r+r1,r),(r,r+r1)]
	
				@test hopp2[:Hopping](args...)[1,1] ≈ HParams.Hopping 
	
			end 
				
		end 
	
	end 
	
end 




#HPerp = TBmodel.BlochHamilt_Perpendicular(latt, hopp2)

#HPar,HPerp = TBmodel.BlochHamilt_ParallPerp(latt, hopp2)

				

#@show Hperp[1]














#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


