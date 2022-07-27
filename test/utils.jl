using myLibs: Utils
using BenchmarkTools






@testset "periodic distance" begin 

	@test Utils.dist_periodic(1,2,1)==0
	@test Utils.dist_periodic(1,1.1,1)≈0.1
	@test Utils.dist_periodic(1.1pi,3pi,2pi)≈0.1pi

	A = rand(2,3) .- 0.5 
	B = rand(2,3) .- 0.5 

	@test Utils.dist_periodic(A,B,0.1,10)≈Utils.dist_periodic(B,A,0.1,10)

	@test Utils.dist_periodic(A[1],B,0.1,10)≈Utils.dist_periodic(B,A[1],0.1,10)
	@test Utils.dist_periodic(A[1],B,0.1,10)[1]≈Utils.dist_periodic(B,A,0.1,10)[1]


end 


@testset "closest data points" begin 

	X = sort(rand(100))

	for i in 5:95 

		i1,i2 = Utils.closest_data_points(X, X[i]/2+X[i+1]/2, 3)

		@test i1 == i-2:i
		@test i2 == i+1:i+3

	end 

end 













error()

















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

result = [1, 2, 4, 6, 7] 

for x in [(5,), (rand(5,3),1), (1:5,), ([1,2,3,4,5],)]

	for y in [(7,), (rand(5,7),2), (1:7,), ([1,2,7],)]

#		println("First: ", typeof.(x)...)
#		println("Second: ", typeof.(y)...)
		
		@assert result==Utils.RescaleInds(x..., y...) 
	
#		println() 
	end 
end 



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


println()
println("----------------")
println("DistributeBallsToBoxes")
println()






for item in Utils.DistributeBallsToBoxes(1,3)

	println(item,"\n")

end 

g14(x) = [x,-x]
f = g14
#Utils.inv_abs 


@time Utils.DistributeBallsToBoxes(2,3,f)
@time Utils.DistributeBallsToBoxes(2,3,f)

@time Utils.DistributeBallsToBoxes(2,3,g14)
@time Utils.DistributeBallsToBoxes(2,3,g14) 

for item in Utils.DistributeBallsToBoxes(2,3,f)

#	println(item," ",sum(abs,item),"\n")
	@assert sum(abs,item) == 2

end 

println("----------------")
println("\nCombsOfVecs\n")

@testset "CombsOfVecs" begin 


vecs= rand(2,5)
coeffs = rand(5,3)


cv = Utils.CombsOfVecs(vecs, coeffs; dim=2)


@test size(cv)==(2,3)



cv2 = Utils.CombsOfVecs(vecs[:,2], coeffs[1:1,:]; dim=2)

@test size(cv2)==(2,3)



cv3 = Utils.CombsOfVecs(vecs, coeffs[:,1:1]; dim=2)


@test size(cv3)==(2,1)




cv4 = Utils.CombsOfVecs(vecs, coeffs[:,1]; dim=2)


@test size(cv4)==(2,) && isapprox(cv4,cv3[:])






end







