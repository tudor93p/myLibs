using myLibs: Utils, ReadWrite

using BenchmarkTools

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




@show Utils.PathConnect(Float64[],100)
@show Utils.PathConnect(Float64[],100,v)



println() 

d = Utils.RecursiveMerge(1:10, dim=2) do i 

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

#	println();	@show n d  

	Utils.nr2string(n,d)# |> println


end 



kPoints = Matrix{Float64}(undef, 0, 1)
NR_KPOINTS = 30

@show Utils.PathConnect(kPoints, NR_KPOINTS, dim=2)




using OrderedCollections: OrderedDict 

Utils.NT(Dict(:a=>2))#|>println

Utils.NT(OrderedDict(:a=>2))#|>println

Utils.NT((a=2,)) #|>println








#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


println()


include("utils_temp1.jl")

@show typeof(X1) 

obs = setdiff(keys(X1),["Energy"]) |> first 
@show obs 
println()

println()

@show size(X1[obs]["A"])
a=rand(3,5);b=rand(3,5);


@show cat(a,b,dims=2) |>size

@show Utils.catNewAx(2,a,b) |> size 
@show Utils.catNewAx(1,a,b) |> size 

@show Utils.catNewAx(2, a,Int.(round.(10b)), rand(3,5,1)) |> size 

println()
Energies = X1["Energy"]*range(1,-1,length=5)

@show Energies


function X(i)

	i==1 && return X1

	X2 = deepcopy(X1)

	X2["Energy"] = Energies[i]

#	for (k,v) in X2[obs]

	return X2 

end 


@show X(1)["Energy"]
@show X(5)["Energy"]


x = Utils.RecursiveMerge(X, 1:5; dim=2)[obs]
	
@show x["A"] |> size


@show X1[obs]["A"] |> size

#
#n = 1000
#
#vecs = map(1:100) do i 
#
#	L = rand(1:10) 
#
#	s = ones(Int,L)
#
#	s[rand(1:L)] = n
#
#	return rand(s...)
#
#end 
#
#Utils.VecsToMat(vecs...; dim=2);
#@time Utils.VecsToMat(vecs...; dim=2);
#@btime Utils.VecsToMat(vecs...; dim=2); 
#



d1 = Dict(1=>2,:a=>:b,3=>rand(2),"s"=>rand(1,1))

d2 = [Dict(1=>3)]

println(Utils.dict_keepkeys(d1,d2...))


@show Utils.uniqlogsp(2,6,6,3; Trunc=true)
#@show Utils.uniqlogsp2(2,6,6,3; Trunc=true)
@show Utils.uniqlinsp(2,6,6,3; Trunc=true)


