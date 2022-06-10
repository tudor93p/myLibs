import myLibs: Utils, Parameters, ComputeTasks 
using OrderedCollections



#===========================================================================#
#
@info "multi task"
#
#---------------------------------------------------------------------------#


module QQ

usedkeys = [:Q1, :Q2]


end 


module RR

usedkeys = [:R12345, :R2]

function adjust_paramcomb(l, P...)

#	@show "here"

	l==2 || return P[l]



	return P[l] 


end 

end 

module SS 

usedkeys = [:S1, :S2]



end 



MDs = map(enumerate([QQ,RR,SS])) do (l,M)

	allparams = Dict() 

	digits = OrderedDict()

	for k in M.usedkeys 

		allparams[k] = 10l .+ sort(Utils.Random_Items(1:10, 0l + rand(2:3)))

		allparams[k] = Utils.flatmap(x->[-x,x,div(x,10),-div(x,10)],allparams[k])|>sort 


		digits[k] = (2,0)

	end 

	return (M,digits,allparams)

end 



function Compute(args...; get_fname::Function, kwargs...)

#	get_fname(args...)() |> println

	A = sum(sum.(values.(args))) .+  rand(2,3)

	sleep(.05)

	return Dict("a"=>A)

end 




println()
													
C = Parameters.Calculation("test_multitask", 
													 Parameters.ParamFlow("Data", MDs...; 
																								adjust=RR.adjust_paramcomb),
													 Compute)  

task0 = ComputeTasks.CompTask(C)


P53 = ComputeTasks.get_rand_paramcomb(task0) 

println("---")
println("---")

task0.get_data(P53...; mute=false)

println("---")

task0.get_data(P53...; mute=true)

println("---")

task0.get_data(P53...)

println("---")
println("---")





for (internal, external) in (
#														([:Q1=>1], Int[]),
#													 ([:Q1=>1], [1,2]),
													 ([:Q1=>1, :R12345=>2], [1,3]),
#													 ([:Q1=>1, :R12345=>2], [2,3]),
#													 ([:Q1=>1, :R12345=>2], [3,4]),
#													 ([:Q1=>1, :Q2=>1,:R12345=>2], [1]),
#													 ([:Q1=>1, :Q2=>1,:R12345=>2], [4]),
													 )

#for e in [(30:40,), external, ((D::AbstractVector{<:Number})->axes(D,1),)] 
for e_ in [(rand(11),rand(7),rand(2))]#, (1,2,3), [(D::AbstractArray{<:Number})->axes(D,i) for i=1:2]]
	
	e = e_[axes(external,1)]


#	second: either vector or function 

#	get_plotaxes(D::AbstractDict)
#
#
#		D -> D["Energy"] 
#
#	end 


@info string("---- ",internal," ------ ", external, " ------ ",typeof.(e))
	


	local out = ComputeTasks.init_multitask(C, internal,
																					[k=>v for (k,v) in zip(external,e)],
																					["Energy", "Atom"]
																					)

	local task, out_dict, construct_Z, = out 
	
	println()
	
	
	for (k,v) in pairs(out_dict)
	
		print(k,"\t",typeof(v),"\t")
	
		v isa AbstractString ? println(v) : println()

	end 

#	println() 


	P = task0.get_plotparams(ComputeTasks.get_rand_paramcomb(task0)...)


#	@show P; 	println()

	P["Q1"]=10

for construct_Z_args in [
												 #(P,),
												 #(P,"some_label"),
												 #(identity, P),
#												 (identity, "a", P),
												 ("a", P),
												 #(identity, P,"some_label"),
												 ]
#	println()
#	@show typeof.(construct_Z_args)

	
	for (k,v) in pairs(construct_Z(construct_Z_args...))

		print(k,"\t",typeof(v),"\t")

#		occursin("line",string(k)) && @show v 

		if v isa AbstractArray 
			
#			print(size(v), "\t", size(v[1]),"\t")
#
#			v[1] isa AbstractArray && print(v[1][1])
#
#			v isa AbstractVector{<:Number} && print(v)

#			println()

		end 

		v isa AbstractString ? println(v) : println()

	end 


	end 

	println()

end end  


#@show Utils.RescaledInds(1:5, 1:10)
#@show Utils.RescaledInds(5, 10) 
#@show Utils.RescaledInds([1,4,5],[1,5]) 








t = ComputeTasks.init_multitask(C, [:Q1=>1], [2=>1], ["Energy"])[1] 


@show t.get_plotparams() 

#t.get_plotparams(Dict()) 
#

t.get_data(Utils.DictFirstVals(t.get_plotparams()),fromPlot=true) 
#t.get_data(Utils.DictFirstVals(t.get_plotparams()),fromPlot=true,mute=false)
#t.get_data(Utils.DictFirstVals(t.get_plotparams()),fromPlot=true,mute=true)


#@show D 

