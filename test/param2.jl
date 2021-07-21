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




allP = Parameters.typical_allparams(((args[:usedkeys], args[:allparams]) for args in [args1,args2])...)


P = Any[1,2]

P[1] = Utils.DictRandVals(allP())

P[2] = Utils.DictRandVals(allP(Dict()))




for (params_digits,p,args) in zip(pds, P,[args1,args2])
	
	@show p

	println("\nParams digits:") 

	foreach(println, params_digits(p))
	
	@show Parameters.tostr(params_digits,p)
	

	println() 




	println()


end

println() 



#ROOT = [ROOT,M1]

get_fname = Parameters.



construct_get_fname(ROOT, ((args[:usedkeys], args[:digits]) for args in (args1, args2))...)


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



a = Parameters.typical_allparams(args1[:usedkeys],args1[:allparams])
b = Parameters.construct_get_fname(ROOT, args1[:usedkeys],args1[:digits])



PF = Parameters.ParamFlow(ROOT, args1[:usedkeys], args1[:digits], args1[:allparams])


@show a() P[1] b(P[1])()
@show PF.get_fname(P[1])()
@show PF.allparams()

println()
a = Parameters.typical_allparams(1, args1[:usedkeys],args1[:allparams])

@show a() 



println()
a = Parameters.typical_allparams(2, args2[:usedkeys], args2[:allparams])
b = Parameters.construct_get_fname(ROOT, args2[:usedkeys], args2[:digits])


@show a(Dict()) 
@show P[2]
@show b(P[2])()


println()

a = Parameters.typical_allparams(((level, args[:usedkeys], args[:allparams]) for (level,args) in enumerate((args1, args2)))...)

@show a() a(Dict()) 

println()


a = Parameters.typical_allparams(((args[:usedkeys], args[:allparams]) for (level,args) in enumerate((args1, args2)))...)
b = Parameters.construct_get_fname(ROOT, ((args[:usedkeys], args[:digits]) for (level,args) in enumerate((args1, args2)))...)

PF = Parameters.ParamFlow(ROOT, ((args[:usedkeys], args[:digits],args[:allparams]) for (level,args) in enumerate((args1, args2)))...)

@show a() a(Dict()) P b(P...)()
@show PF.get_fname(P...)()
@show PF.allparams(P[1:end-1]...)

println()

a = Parameters.typical_allparams(2, ((level, args[:usedkeys], args[:allparams]) for (level,args) in enumerate((args1, args2)))...)



@show a() a(Dict()) 

println()
a = Parameters.typical_allparams(((level, args[:usedkeys], args[:allparams]) for (level,args) in enumerate((args1, args2)))...)


@show a() a(Dict()) 

println()

@show Parameters.good_methods(1, identity)
@show Parameters.good_methods(nothing, identity)
@show Parameters.good_methods(PF, merge)
@show Parameters.good_methods(PF, identity) 








#construct_get_fname(path, (uk1, d1), (uk2, d2), ...)

#typical_allparams((uk1, ap1), (uk2, ap2), ...)


#ParamFlow(NrParamSets, path, (uk1, d1, ap1), (uk2, d2, ap2), ...)




#uk1,d1,ap1 
#uk1,d1,a1 
#pd1,uk1,ap1 
#pd1,a1 
#
#
#(uk1,d1) or pd1 
#(uk2,d2) or pd2
#(uk1,ap1) or a1 
#(uk2,ap2) or a2 






















