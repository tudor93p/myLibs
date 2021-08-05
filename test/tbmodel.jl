using myLibs: TBmodel, Operators, Lattices 
using LinearAlgebra 


R = rand(2,3);
R= hcat(R,R.+[0,1],R.+[1,0])#,R.+[0,-1],R.+[-1,0])

hopp(ri,rj) = isapprox(norm(ri-rj),0)*2.0+isapprox(norm(ri-rj),1)*1.0




#println(TBmodel.HoppingMatrix(R; Hopping=hopp))

TBL = ([0,1], [0,[2,0]], R)

inter, (ms,Rms,Tms) = TBmodel.Compute_Hopping_Matrices(TBL,Hopping=hopp)


@show size(inter)

@show ms 
@show Rms 

@show size.(Tms)



H = TBmodel.Bloch_Hamilt(TBL; Hopping=hopp)


@show H(1) - H([2]) + H([3,4]) - H([5,6,7])



println()
println()
println()

F = Operators.Position_Expectation(1, R, nr_orb=1)#; kwargs...)#, fpos=df[2])


E =  LinearAlgebra.eigen(Matrix(H(1)))

#@show propertynames(E)


P = E.vectors;


@assert isapprox(H(1)*P[:,1:1] ,E.values[1]*P[:,1:1])


@show size(F(P[:,1:3]))
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
@show Lattices.VecDim(L)
@show Lattices.LattDims(L)

println() 

tbl = Lattices.NearbyUCs(L)

println.(tbl)


intra,inter = TBmodel.Compute_Hopping_Matrices(tbl; Hopping=hopp)

println()


@show intra 
println("\ninter:")
!isnothing(inter) && println.(inter)

println()



H = TBmodel.Bloch_Hamilt(tbl; Hopping=hopp)


@show H(rand(Lattices.LattDim(L)))


println()

Hpar = TBmodel.BlochHamilt_Parallel(L, Dict(:Hopping=>hopp))


@show Hpar(rand(Lattices.LattDim(L)-1))

println()



Hperp = TBmodel.BlochHamilt_Perpendicular(L, Dict(:Hopping=>hopp))



@show Hperp[1]




HPar,Hperp = TBmodel.BlochHamilt_ParallPerp(L, Dict(:Hopping=>hopp))




@show Hperp[1]


@show Hpar(rand(Lattices.LattDim(L)-1))














