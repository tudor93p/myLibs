using myLibs: Parameters ,Utils
using OrderedCollections: OrderedDict 
using Combinatorics: powerset

input_dict = Dict(

		:allparams => Dict(:X=>[1],:Y=>[3,4],:T=>[10,20,30]),

#		:digits => OrderedDict(:X=>(1,0), :Y=>(1,1), :Z=>(2,2)),
		:digits => (X=(1,0), Y=(1,1), Z=(2,2)),

		)





module M 

	import ..input_dict,..Parameters
	#export usedkeys
	usedkeys(args...) = [:X]
	a=3 
	
	NrParamSets = 2 
	
	digits = (X=(3,3),) 
	
	path = "abc" 
	
	params_digits = Parameters.typical_params_digits(usedkeys, input_dict[:digits])
	
	
	
	

end 

function usedkeys()
	[:X,:Z,:zze]

end 


usedkeys(P::AbstractDict) = setdiff(usedkeys(),[:zze])

NrParamSets = 1

usedkeys_symb = [:Y,:T]

usedkeys22() = [:X]

usedkeys22(P::AbstractDict) = usedkeys22()



Parameters.typical_allparams(usedkeys, input_dict[:allparams]) |> println
Parameters.typical_allparams(usedkeys_symb, input_dict[:allparams]) |> println

Parameters.typical_allparams(usedkeys22, input_dict[:allparams]) |> println



#@show Parameters.get_usedkeys(usedkeys)
#@show Parameters.get_usedkeys(usedkeys22)
#@show Parameters.get_usedkeys(usedkeys_symb)

P = Dict(:X=>44,:Z=>0.5,:Y=>-3,:T=>-2.5) 


@show Parameters.get_usedkeys(usedkeys, P)
@show Parameters.get_usedkeys(usedkeys22, P)
@show Parameters.get_usedkeys(usedkeys_symb, P)



@show Parameters.params_digits_(P, usedkeys, input_dict[:digits])
#@show Parameters.params_digits_((P,P), usedkeys22, input_dict[:digits])
@show Parameters.params_digits_(Dict(), usedkeys_symb, input_dict[:digits])


println();















ROOT = "Data"

#jFN = Parameters.FileName("Data", usedkeys, identity, identity) 


#@show Parameters.fname(ROOT)("some_file")

#@show FN


pd = Parameters.typical_params_digits(usedkeys, input_dict[:digits])

@show pd(P)

println()

@show Parameters.params_to_string(pd(P))

@show Parameters.prefix(ROOT)
@show Parameters.prefix(M)
@show Parameters.prefix(ROOT,M)
@show Parameters.prefix(M,ROOT)



#
#@show Parameters.fname(ROOT, pd)("some_file")
#@show Parameters.fname(ROOT, pd, P)("some_file")
#@show Parameters.fname(ROOT, pd, P, P)("some_file")
#
#
#@show Parameters.fname(M, pd)("some_file") 
#@show Parameters.fname(M, pd, P)("some_file") 
#@show Parameters.fname(ROOT, M, pd, P)("some_file") 
#@show Parameters.fname(ROOT, M, pd, P)("some_file") 
@show Parameters.prefix("", M,ROOT)#, pd, P)("some_file") 
#
#
##ROOT, M, pd
#
#
#get_fname(args...) = Parameters.fname(ROOT, M, pd, args...) 
#
#@show get_fname(P,P,P)("nnewfile")
#

println()


pd = Parameters.typical_params_digits(usedkeys, input_dict[:digits])

#for pd_ in [pd, input_dict[:digits]] 

for path in [ROOT, (ROOT,M)]

#	for args in [ Any[usedkeys, path, pd], Any[usedkeys, input_dict[:digits], path]]

for args in [Any[usedkeys, input_dict[:digits]], 
						 Any[pd],
						 Any[M],
						 Any[M, input_dict[:digits]]
						 ]



#	for inds in powerset(1:length(args)) 
	
		args_ = copy(args) 
	
#		for i in inds 
	
#			args_[i] = M
	
#		end 
	
#		@info "args_"

#		foreach(println,args_)
#		println()
#		local FN = Parameters.FilenameGenerator(args_...)
	
#		@show FN.get_fname(P)("f")

		
	@show Parameters.construct_get_fname(path, args_...)(P)("f")



#		local PF = Parameters.ParamFlow(identity, args_...)
#
#		PF.get_fname(P)("f")==FN.get_fname(P)("f") || error()
#
#
#		local PF = Parameters.ParamFlow(input_dict, args_...)
#
#		@show PF.allparams() 
#
#		println(Parameters.ParamFlow((M,input_dict[:allparams]), 
#																 args_...).allparams())
		println()



	end 

end 
#end 

for a1 in [usedkeys,M]

	@show Parameters.typical_allparams(a1, input_dict[:allparams])
end 


#usedkeys, allp 



#FN = Parameters.FilenameGenerator(ROOT, pd)


#@show FN.get_fname(P)("f") 

println() 










input_dict = Dict(:allparams=>(
										length = [10,20],
									 	width = [7],
										Barrier_height = [0,0.5],
										SCpx_magnitude = [0.6],
										),

									 :digits=>(
											length = (3, 0),

											Barrier_height = (1,3),
										)	
									
									)

usedkeys2 = [:length, :width, :SCpx_magnitude]

#FN = Parameters.FilenameGenerator(usedkeys2, input_dict[:digits], "Data")


#@show FN.usedkeys 


#PF = Parameters.ParamFlow(input_dict[:allparams], usedkeys2, input_dict[:digits], "Data")
#
#
#@show PF.allparams()
#
#PF = Parameters.ParamFlow(identity, usedkeys2, input_dict[:digits], "Data")
#
#@show PF.allparams(0)
#
#
#PF = Parameters.ParamFlow(1, input_dict[:allparams], usedkeys2, 
#													input_dict[:digits], "Data")
#
#
#@show PF.NrParamSets
#
#
#
#@show PF.allparams()
#
#
# 
#
#


@show Parameters.union_usedkeys(usedkeys2, usedkeys)(P)  
@show Parameters.union_usedkeys(usedkeys2, usedkeys2)

println() 

#@show Parameters.union_usedkeys(usedkeys, PF)(P)
#@show Parameters.union_usedkeys(usedkeys, PF)()




println()

#for item  in Parameters.combine_params_digits(PF, PF.params_digits, (P, Symbol[:x],(x=(1,2),)))
#
#	println(item)
#
#end 
#
#println()
#
#@show Parameters.combine_params_digits(PF, PF, PF)(P)
#
#
println()




#
#module M3
#
#
#Read(P;some_kwarg=2,kw...) = ("read",some_kwarg)
#
#Compute(P;some_kwarg=1,kw...) = ("compute", some_kwarg)
#
#FoundFiles(P;kwa...) = rand(Bool)
#
#NrParamSets = 1
#	
#apnt =(
#										length = [10,20],
#									 	width = [7],
#										Barrier_height = [0,0.5],
#										SCpx_magnitude = [0.6],
#										)  
#
#allparams() = Dict(k=>apnt[k] for k in keys(apnt))
#
#end
#
#
# 
#
#@show Parameters.Calculation("calc1", PF, M3).Compute(P)
#@show Parameters.Calculation(PF, M3; some_kwarg=3).Compute(P)
#@show Parameters.Calculation(PF, M3).Compute(P, some_kwarg=4)
#
#
#
#println()
#
#
#C = Parameters.Calculation("calc2", PF, M3) 
#
#@show C.name
#@show Parameters.Calculation(PF, M3).name
#
#@show Parameters.get_NrPSets(C)
#@show Parameters.get_NrPSets(PF)
#
#@show Parameters.get_allP(C)
#@show Parameters.get_allP(PF)
#
#
#
#println() 
#
#foreach(println, Parameters.get_paramcombs(C))
#
#
#println() 
#
#foreach(println, Parameters.get_paramcombs(PF))
#
#
#
println()




println() 
















 















 















 















 















 















 















 















 















 















 















 















 















 















 















 















 















 















 















 















 















 















 















 















 















 















 















 















 















 















 















 















 















 















 















 















 















 















 















 















 















 















 















 















 















 















 















 















 















 















 















 















 















 















 















 















 















 















 















 















 















 















 















 















 















 















 















 















 















 















 















 















 















 















 















 















 















 















 















 















 















 















 















 















 















 















 















 















 















 















 













 















 















 















 















 




















println()
