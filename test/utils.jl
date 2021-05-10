

using myLibs: Utils 


P = rand(10,7) 


p1 = Utils.PathConnect(P, 100)[1]


@show size(p1)


p2 = Utils.PathConnect(transpose(P), 100,dim=2)[1]


@show size(p2)




@show isapprox(p1,transpose(p2))


#p1 .- transpose(p2) 

#@show P[3,:] p1[end,:] p2[:,end]














































nothing
