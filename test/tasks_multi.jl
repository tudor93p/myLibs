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

		allparams[k] = 10l .+ sort(Utils.Random_Items(1:10, rand(l:l+5)))

		digits[k] = (2,0)

	end 

	return (M,digits,allparams)

end 



function Compute(args...; get_fname::Function)

	@show get_fname(args...)()

	return sum(sum.(values.(args))) .+  rand(20)

end 




													
C = Parameters.Calculation("test_multitask", 
													 Parameters.ParamFlow("Data", MDs...),
													 Compute) 


task, out_dict, construct_Z, = ComputeTasks.init_multitask(
																	C,
																	[:Q1=>1], # internal keys 
																	[2=>1:20]; # external param
																	)












