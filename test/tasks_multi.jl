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

function adjust_paramcomb(l, P...)

	@show "here"

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

		allparams[k] = 10l .+ sort(Utils.Random_Items(1:10, 0l + rand(1:2)))

		digits[k] = (2,0)

	end 

	return (M,digits,allparams)

end 



function Compute(args...; get_fname::Function, kwargs...)

#	get_fname(args...)() |> println

return sum(sum.(values.(args))) .+  rand(5)

end 




println()
													
C = Parameters.Calculation("test_multitask", 
													 Parameters.ParamFlow("Data", MDs...; 
																								adjust=RR.adjust_paramcomb),
													 Compute) 

for (internal, external) in (
#													 ([:Q1=>1], [1,2]),
#													 ([:Q1=>1, :R1=>2], [1,3]),
#													 ([:Q1=>1, :R1=>2], [2,3]),
#													 ([:Q1=>1, :R1=>2], [3,4]),
#													 ([:Q1=>1, :Q2=>1,:R1=>2], [1]),
#													 ([:Q1=>1, :Q2=>1,:R1=>2], [4]),
													 )

#for e in [(30:40,), external, ((D::AbstractVector{<:Number})->axes(D,1),)] 
for e in [(rand(11),rand(7)), (1,2), [(D::AbstractArray{<:Number})->axes(D,i) for i=1:2]]
	



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


	P = task.get_plotparams(rand(task.get_paramcombs())...)

#	@show P; 	println()


for construct_Z_args in [
												 (P,),
												 (P,"some_label"),
												 (identity, P),
												 (identity, P,"some_label"),
												 ]
#	println()
#	@show typeof.(construct_Z_args)

	for (k,v) in pairs(construct_Z(construct_Z_args...))

		print(k,"\t",typeof(v),"\t")

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


for P in t.get_paramcombs()

	println(P,"\n")
end


t.get_plotparams() |> println
