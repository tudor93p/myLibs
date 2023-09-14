import myLibs: Utils, Geometry

using OrderedCollections:OrderedDict 

include("../src/TightBinding/Lattices.jl")

include("../src/TightBinding/Lattices_VecMatDict.jl")




A = rand(2,3)

A[:,1] .= 1 
A[:,2] .= 2
A[:,3] .= 3



@show NrVecs(A)
@show NrVecs(A,nothing) 
@show NrVecs(A,3) 

@show IndsVecs(A)
@show IndsVecs(A,1)
@show IndsVecs(A,1:2)

@show Vecs(A)
@show Vecs(A,2)
@show Vecs(A,[1,3]) |> size  
@show Vecs(A,[1 3])  |> size 


@show VecLen(A)
@show VecAx(A) 

@show EmptyVec()
@show EmptyVec(3) 
@show EmptyVec(A) 


@show Cat(A,EmptyVec(A)) 

@show ZeroVecs(2,3) == zero(A) 



l2 = to_myODict(["A"=>rand(2,1),rand(2,5),Dict(2=>rand(2,3)),(["A","B"],rand(2,2)),rand(2,4)],2) 



#@show l2 
