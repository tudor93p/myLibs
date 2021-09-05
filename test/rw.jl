import myLibs:ReadWrite, Utils





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




