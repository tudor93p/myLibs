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

usedkeys = [:R1, :R2]


end 

module SS 

usedkeys = [:S1, :S2]



end 



MDs = map(enumerate([QQ,RR,SS])) do (l,M)

	allparams = Dict() 

	digits = OrderedDict()

	for k in M.usedkeys 

		allparams[k] = 10l .+ sort(Utils.Random_Items(1:10, l + rand(1:5)))

		digits[k] = (2,0)

	end 

	return (M,digits,allparams)

end 



function Compute(args...; get_fname::Function, kwargs...)

#	get_fname(args...)() |> println

	return sum(sum.(values.(args))) .+  rand(10)

end 




println()
													
C = Parameters.Calculation("test_multitask", 
													 Parameters.ParamFlow("Data", MDs...),
													 Compute) 

for (internal, external) in (
													 ([:Q1=>1], [1]),
													 ([:Q1=>1, :R1=>2], [2]),
#													 ([:Q1=>1, :R1=>2], [3]),
#													 ([:Q1=>1, :Q2=>1,:R1=>2], [1]),
#													 ([:Q1=>1, :Q2=>1,:R1=>2], [4]),
													 )

for e in [(30:40,), ((D::AbstractDict)->40:50,)] 
	



#	second: either vector or function 

#	get_plotaxes(D::AbstractDict)
#
#
#		D -> D["Energy"] 
#
#	end 


	@info string("---- ",internal," ------ ", external)
	
	local out = ComputeTasks.init_multitask(C, internal,
																					[k=>v for (k,v) in zip(external,e)]
																					)

	local task, out_dict, construct_Z, = out 
	
	println()
	
	
	for (k,v) in pairs(out_dict)
	
		println(k,"\t",v)
	
	end 

#	println() 


	P = task.get_plotparams(rand(task.get_paramcombs())...)

#	@show P; 	println()


for construct_Z_args in [
												 (P,),
#												 (P,"some_label"),
#												 (identity, P),
#												 (identity, P,"some_label"),
												 ]
#	println()
#	@show typeof.(construct_Z_args)

	for (k,v) in pairs(construct_Z(construct_Z_args...))

		print(k,"\t",typeof(v),"\t")

		if v isa AbstractArray
			
			print(size(v), "\t", size(v[1]),"\t")

			v[1] isa AbstractMatrix && print(v[1][1])

			println()

		end 

		v isa AbstractString && println(v)

	end 


	end 

	println()

end end  


#@show Utils.RescaledInds(1:5, 1:10)
#@show Utils.RescaledInds(5, 10) 
#@show Utils.RescaledInds([1,4,5],[1,5]) 









