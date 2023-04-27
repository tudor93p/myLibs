import myLibs: TBmodel, Operators, Lattices, H_Superconductor, ArrayOps
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

