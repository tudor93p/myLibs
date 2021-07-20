import myLibs: Parameters, Utils 



ROOT = "Data"

#===========================================================================#
#
# Level 1
#
#---------------------------------------------------------------------------#

module M1 

NrParamSets = 1

export lev1_input 

lev1_input = Dict(


			:allparams => Dict(:X1 => [1,2],
												 :Y1 => [-0.5,0.9,1.3],
												 :Z1 => [-10,20,35],
												 :T1 => ("k",8),
												 ),

			:digits => (X1 = (1,1), Y1=(1,2), T1=())

			)


lev1_usedkeys = [:X1,:T1,:Z1] 

end 



















#===========================================================================#
#
# Levels 1 and 2
#
#---------------------------------------------------------------------------#

module M2 

using ..M1: lev1_input


NrParamSets = 2 

lev2_input = Dict(


			:allparams => Dict(:X2 => [-5,100],
												 :Y2 => [-0.05,0.009,1.003],
												 :Z2 => [-10,20,35],
												 ),

			:digits => (Z2 = (3,0), Y2=(1,3))

			)


lev2_usedkeys = [:Z2,:Y2,:X2]

end 

#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


args1 = merge(M1.lev1_input,
							Dict(:usedkeys => M1.lev1_usedkeys, :path => ROOT,)
							)

args2 = merge(M2.lev2_input,
							Dict(:usedkeys => M2.lev2_usedkeys, :path => ROOT,)
							)


pds = Parameters.typical_params_digits(((args[:usedkeys], args[:digits]) for args in (args1, args2))...)


#allP = [Parameters.typical_allparams(args[:usedkeys], args[:allparams]) for args in [args1,args2]]


allP = Parameters.typical_allparams(((args[:usedkeys], args[:allparams]) for args in [args1,args2])...)



P = Utils.DictRandVals.(allP)



for (params_digits,p,allp,args) in zip(pds, P,allP,[args1,args2])
	@show allp
	
	@show p

	println("\nParams digits:") 

	foreach(println, params_digits(p))
	
	@show Parameters.tostr(params_digits,p)
	

	println() 




	println()


end

println() 



#ROOT = [ROOT,M1]

get_fname = Parameters.construct_get_fname(ROOT, ((args[:usedkeys], args[:digits]) for args in (args1, args2))...)


fname = get_fname(P...)


@show fname("test")

println()
println()

get_fname = Parameters.construct_get_fname(ROOT, args1[:usedkeys], args1[:digits])



fname = get_fname(P[1]) 

@show fname("test")


println()
println()



#$get_fname = Parameters.construct_get_fname(ROOT, M1, M1.lev1_input)



#@show get_fname(P[1])("test")





allP = Parameters.typical_allparams(((args[:usedkeys], args[:allparams]) for args in (args1, args2))...)


@show allP





#construct_get_fname(path, (uk1, d1), (uk2, d2), ...)

#typical_allparams((uk1, ap1), (uk2, ap2), ...)











