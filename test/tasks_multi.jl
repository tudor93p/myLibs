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

		allparams[k] = 10l .+ sort(Utils.Random_Items(1:20, l + rand(10:10)))

		digits[k] = (2,0)

	end 

	return (M,digits,allparams)

end 



function Compute(args...; get_fname::Function, kwargs...)

#	get_fname(args...)() |> println

	return sum(sum.(values.(args))) .+  rand(100)

end 




println()
													
C = Parameters.Calculation("test_multitask", 
													 Parameters.ParamFlow("Data", MDs...),
													 Compute) 

for (internal, external) in (
#													 ([:Q1=>1], [1]),
#													 ([:Q1=>1, :R1=>2], [2]),
#													 ([:Q1=>1, :R1=>2], [3]),
													 ([:Q1=>1, :Q2=>1,:R1=>2], [1]),
#													 ([:Q1=>1, :Q2=>1,:R1=>2], [4]),
													 )

	@info string("---- ",internal," ------ ", external)
	
	local out = ComputeTasks.init_multitask(C, internal, external)
	
	local task, out_dict, construct_Z, = out 
	
	println()
	
	
	for (k,v) in pairs(out_dict)
	
		println(k,"\t",v)
	
	end 

	println()
end 


#@show Utils.RescaledInds(1:5, 1:10)
#@show Utils.RescaledInds(5, 10) 
#@show Utils.RescaledInds([1,4,5],[1,5]) 









