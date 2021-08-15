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



function Compute(args...; get_fname::Function)

	@show get_fname(args...)()

	return sum(sum.(values.(args))) .+  rand(20)

end 




println()
													
C = Parameters.Calculation("test_multitask", 
													 Parameters.ParamFlow("Data", MDs...),
													 Compute) 


out = ComputeTasks.init_multitask(C,
																	[:Q1=>1], #internal key => x
#																	[:Q1=>1, :R1=>2], # x and z 
																	[2]; # external param =>y
																	)


task, out_dict, construct_Z, = out 

println()


for (k,v) in pairs(out_dict)

	println(k,"\t",v)

end 




#@show Utils.RescaledInds(1:5, 1:10)
#@show Utils.RescaledInds(5, 10) 
#@show Utils.RescaledInds([1,4,5],[1,5]) 









