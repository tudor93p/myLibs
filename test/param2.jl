import myLibs: Parameters, Utils 
using myLibs.Parameters


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


									:allparams => Dict{Symbol,Any}(:X1 => [1,5],
																		 :Y1 => [-0.5,0.9],#,1.3],
																		:Z1 => [-10,20],#,35],
												 :T1 => ("k",8),
												 ),

									:digits => (X1 = (1,1), Y1=(1,2), T1=(),Q1=(1,1))

			)


lev1_usedkeys = [:X1,:T1,:Z1,:Q1]
usedkeys=lev1_usedkeys



function adjust_paramcomb(l::Int, P...)

	l==1 || return P[l]

	p = copy(P[1])

	isa(p[:X1], Number) && return setindex!(p, p[:X1]/2, :Q1)
	
	return delete!(p, :Q1)

end






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


									:allparams => Dict{Symbol,Any}(:X2 => [-5,100],
																		 :Y2 => [-0.05],#,1.003],
												 :Z2 => [20],#,35],
												 ),

									:digits => (Z2 = (3,0), Y2=(1,3), Label=())

			)


lev2_usedkeys = [:Z2,:Y2,:X2,:Label]

usedkeys=lev2_usedkeys 



function adjust_paramcomb(l, P...) 
	
	isa(P[1][:X1],Number) && @assert haskey(P[1], :Q1)

	return P[l]
end 

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

#if false   ##############3

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


fname72 = get_fname(P...)


@show fname72("test")

println()
println()

get_fname = Parameters.construct_get_fname(ROOT, args1[:usedkeys], args1[:digits])



fname72 = get_fname(P[1]) 

@show fname72("test")


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





#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

println()
println()


a = Parameters.typical_allparams(args1[:usedkeys], [args1[:allparams] for i=1:2])

b = Parameters.construct_get_fname(ROOT, args1[:usedkeys], [args1[:digits] for i=1:2])


@show a()																 

p = Utils.DictRandVals(a())

@show p


@show b(p)()



println()


a = Parameters.typical_allparams(((args[:usedkeys], [args[:allparams] for i=1:2]) for args in (args1,args2))...)
b = Parameters.construct_get_fname(ROOT, ((args[:usedkeys], [args[:digits] for i=1:2]) for args in (args1,args2))...)
																 
																 
@show a() a(Dict())



P = Utils.DictRandVals([a(), a(Dict())])

@show P  length.(P)


@show b(P...)()
println()





a = Parameters.typical_allparams((args1[:usedkeys], args1[:allparams]), (args2[:usedkeys], [args2[:allparams] for i=1:3]))
b = Parameters.construct_get_fname(ROOT, (args1[:usedkeys], args1[:digits]), (args2[:usedkeys], [args2[:digits] for i=1:3]))

function cond(level::Int, P...)::Bool

	level==2 || return true 

	v1 = P[1][:Z1]
	
	v2, v3 = getindex.(P[2], :X2)

	any(v->v isa AbstractVector, (v1,v2,v3)) && return true 

	return v1*v2*v3<0

end  

function cond1(level::Int, p) 
	
	level==1 && return p[:X1]>3

	return true 

end 


function adjust(level::Int, P...)

#	p = copy(P[level])

#	if level==1 
#
#		newv(x::Number) = x < 0 ? abs(x)+5 : x
#
#		newv(a::AbstractVector) = newv.(a) 
#
#		p[:Z1] = newv(p[:Z1])
#
#	end 

	if level==2 

#		@show getindex.(P[level],:X2)

#@show P[1] P[2] length(P)

	@assert Utils.isList(P[2],UODict) P[2]

		newv(x::Number) = x<0 ? abs(x) + abs(P[1][:X1]) + rand() : x

		newv(x::AbstractVector) = newv.(x)

#		println("I'm adding :X7")
		

		for (i,p_) in enumerate(P[level])
	
			P[level][i][:X2] = newv(p_[:X2])
			
#			P[level][i][:X7]=7
		end 

#		@show getindex.(P[level],:X2)
		
	end 

	return P[level]

end 



@show adjust 


c = Parameters.ParamFlow(ROOT, 
												 (args1[:usedkeys], args1[:digits],args1[:allparams]), 
												 (args2[:usedkeys], [args2[:digits] for i=1:3], [merge(args2[:allparams],Dict(:Label=>string(l))) for l='A':'C']),
												 constraint=cond1, adjust=adjust
												 )

																 
@show a() a(Dict())
@show c.allparams() c.allparams(Dict())


println()

@info "All combinations (a few shown randomly)"



for comb in (Parameters.get_paramcombs(c) |> C -> Utils.Random_Items(C, min(3,length(C))))

	foreach(println, comb)

	println() 

end 

println() 

@info "Choose one comb"

P = rand(Parameters.get_paramcombs(c))


println("####################################")

println("P");println.(P)

println()


@show length.(Utils.flat.(P))


@show b(P...)()
@show c.get_fname(P...)()


println() 


@show c.NrParamSets

@show Parameters.convertParams_toPlot(c)



plot_P  = Parameters.convertParams_toPlot(c, P)


println() 

@show plot_P



P_  = Parameters.convertParams_fromPlot(c, plot_P)


@show P_  

function test(p::AbstractDict, p_::AbstractDict) 
	
	Set(keys(p))==Set(keys(p_)) || return false 

	for (k,v) in pairs(p)

		p_[k]==v || return false 

	end 

	return true 

end 


test(p::Utils.List, p_::Utils.List) = all(test.(p,p_))
#test((x,y)) = test(x,y)



@show  [test(p...) for p in zip(P,P_)]  #|> all



#end # first 'if'


println("############################")



c = Parameters.ParamFlow(ROOT, 
												 (args1[:usedkeys], args1[:digits], args1[:allparams]), 
												 (args2[:usedkeys], [args2[:digits] for i=1:2], [merge(args2[:allparams],Dict(:Label=>string(l))) for l='A':'B']),
												 adjust=(M1.adjust_paramcomb,M2.adjust_paramcomb)
												 )





@show c.allparams()
@show c.allparams(Dict())



Parameters.convertParams_toPlot(c) |>println



#error() 

pickfirst(l,P...)=Utils.DictFirstVals(P[l])

for item in Parameters.get_paramcombs(c; repl=pickfirst)
																			
	println("\n ---------- Original P -----------\n")
	println.(item)
	
	println("\n ---------- Plot P -----------\n")

	pv = Parameters.convertParams_toPlot(c, item)

	println(pv)

	println("\n ---------- Recovering P -----------\n")


	item2 = Parameters.convertParams_fromPlot(c, pv)

	println.(item2)

	println()

@show  [test(p...) for p in zip(item,item2)]  #|> all
	break

end 









































