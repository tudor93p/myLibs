

using myLibs: Utils, ReadWrite


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
			 "x"=>rand(1:10)+i,
			 "y"=>rand(1:10,3).+i,
			 "z"=>Dict("z1"=>rand(1:10),"z2"=>rand(1:10,2)), 	
#			 "t"=>OrderedDict("t1"=>Dict(1=>rand()), "t2"=>3,"t3"=>rand(5))
			 )

end 


@show d["x"] |> size  
println()
@show d["y"] |> size 
println()
@show d["z"]["z1"] |> size 
@show d["z"]["z2"] |> size 
#@show d["t"] |> typeof
#@show d["t"]["t1"] |> typeof
#

println()
println()

FSM = "dat"

fname = x->string("test/savefile/",x)

Write!,outdict = ReadWrite.Write_PhysObs(fname, FSM)

for k in keys(d) 

@show k 

#foreach(println∘typeof,values(d[k]))





Write!(k, d[k], outdict)

z = ReadWrite.Read_PhysObs(fname, k, FSM)[k] 

println(d[k])
println(z) 


println()

end 




@show Utils.IdentifySectors([1,1,1,3,4,5,5,1]) 

@show Utils.IdentifyRanges([1,1,2,3,4,5,5,1])

@show Utils.IdentifyRanges([1,1,3,5,7,5,5,1,2,3])

@show Utils.IdentifyRanges(Utils.Random_Items(1:10))



for n in ["abc",1,10.5], d in [(),1,3,(1,2),(2,)]

	@show n d  

	println(Utils.nr2string(n,d))

	println()

end 



kPoints = Matrix{Float64}(undef, 0, 1)
NR_KPOINTS = 30

@show Utils.PathConnect(kPoints, NR_KPOINTS, dim=2)




using OrderedCollections: OrderedDict 

Utils.NT(Dict(:a=>2))|>println

Utils.NT(OrderedDict(:a=>2))|>println

Utils.NT((a=2,)) |>println









nothing
