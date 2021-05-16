

using myLibs: Utils 


P = rand(10,7) 


p1 = Utils.PathConnect(P, 100, dim=1)[1]


@show size(p1)


p2 = Utils.PathConnect(transpose(P), 100,dim=2)[1]


@show size(p2)




@show isapprox(p1,transpose(p2))


#p1 .- transpose(p2) 

#@show P[3,:] p1[end,:] p2[:,end]




n = 123 
v = [1,0,3,1]

d =  Utils.PropDistributeBallsToBoxes(n, v)

@show d sum(d) n Utils.LA.normalize(d,1)*sum(v)


v = [1,1,3,1]

println.(Utils.PathConnect([-1.2,-0.6,0,0.6,1.2],10))

println() 

println.(Utils.PathConnect([-1.2,-0.6,0,0.6,1.2],10,v))




@show Utils.PathConnect([],100)
@show Utils.PathConnect([],100,v)



println() 

d = Utils.RecursiveMerge(1:10) do i 

		Dict(	
			 "x"=>rand()+i,
			 "y"=>rand(3).+i,
			 "z"=>Dict("z1"=>rand(),"z2"=>rand(2)), 	
			 "t"=>OrderedDict("t1"=>Dict(1=>rand()), "t2"=>3)
			 )

end 


@show d["x"] |> size  
println()
@show d["y"] |> size 
println()
@show d["z"]["z1"] |> size 
@show d["z"]["z2"] |> size 
@show d["t"] |> typeof
@show d["t"]["t1"] |> typeof
























nothing
