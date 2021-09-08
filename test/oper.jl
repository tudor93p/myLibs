import myLibs:Operators ,Utils
import LinearAlgebra 


nr_at = 10

nr_orb = 3 
nr_wf = 17 

size_H = nr_at*nr_orb 


E = sort(rand(nr_wf))*2 .- 1 

dim = 2

dim2 = [2,1][dim]

atoms = rand([nr_at,2][[dim,dim2]]...)

shape = [0,0]

shape[dim]=nr_wf
shape[dim2]=size_H


P = mapslices(LinearAlgebra.normalize,rand(ComplexF64,shape...),dims=dim2) 


for psi in eachslice(P,dims=dim)

	@assert isapprox(1,LinearAlgebra.norm(psi))

end  

kwargs = (dim=dim,nr_at=nr_at,nr_orb=nr_orb,size_H=size_H)


for Op in [#Operators.Position_Expectation(1, atoms; kwargs...),
					 #Operators.Operator([-1,0,1]; kwargs...),
					 #Operators.IPR(;kwargs...),
					 #Operators.LDOS([1]; kwargs...),
					 ]


	@show Op


end 


for n in 0:3

	@show n

	ks = Utils.Random_Items([:nr_orb,:nr_at,:size_H],n)

	@show ks 


	kw = Utils.dict_diff(kwargs,ks)

	@show kw 

	@show Operators.get_nratorbH(;kw...) 

#	if n>0
#
#		k = rand(ks) 
#	
#		@show k 
#		@show Operators.get_nratorbH(kwargs[k];kw...) 
#	
#	end 


i = rand(2:length(kwargs))

	@show i 

	k = keys(kwargs)[i]

	@show k
	v = kwargs[k]
@show v 

	A = rand(v)


					 

	@show size(A)

	@show Operators.Operator(A; kwargs...)

	println()
	println()
end  




for k in [:nr_orb,:nr_at,:size_H]

	nrs,enough_data = Operators.get_nratorbH(;kwargs...)

	@show k Operators.which_diag(rand(kwargs[k]),nrs)

end 

println()
println()

function po(Op)

	for k in propertynames(Op)

		k in [:data, :diag, :acts_on_atoms, :acts_on_orbs, :inds] || continue 

		v = getproperty(Op, k)


		print(k," => ") 
		
		if k == :data 
			
			println(typeof(v)," ",size(v))

		else 

			println(v) 

		end 
	
	end
	
	println()

end 


#
#println("\n------------ LDOS ---------------\n")
#
#Op = Operators.LDOS(;kwargs...)
#
#po(Op)
#
#
#@time v = Operators.DEF(P, Op)
#@time v = Operators.DEF(P, Op)
#
#@show size(v)
#
#
#
#println("\n------------- Position 1 --------------\n")
#
#
#Op = Operators.Position(1, atoms; kwargs...)
#
#po(Op)
#
#@time v = Operators.DEF(P, Op)
#@time v = Operators.DEF(P, Op)
#
#@show size(v)
#
#@show v[5] v[9]
#
#println("\n------------ Position 1 no sum_atoms ---------------\n")
#
#Op = Operators.Position(1, atoms; sum_atoms=false, kwargs...)
#
#po(Op)
#@time v = Operators.DEF(P, Op)
#@time v = Operators.DEF(P, Op)
#
#@show size(v)
#
#
#
#@show sum(v[:,5]) sum(v[:,9])
#
#
#
#println("\n------------ Position 1 no sum_orbs ---------------\n")
#
#Op = Operators.Position(1, atoms; sum_orbitals=false, kwargs...)
#
#po(Op)
#@time v = Operators.DEF(P, Op)
#@time v = Operators.DEF(P, Op)
#
#@show size(v)
#
#
#
#@show sum(v[:,5]) sum(v[:,9])


for n in [1,nr_at,nr_orb,size_H] 
	
for N in [(n,),]#(n,n)]

#for n in [1,nr_at], N in [(n,),(n,n)]
		
	R = round.(rand(N...), digits=2) 

	X = [] 

	for (sum_atoms,sum_orbitals) in [(true,true),
																	 (true,false),
																	 (false,true),
																	 (false,false)
																	 ]
		println("\n---------------------------")

		@show N sum_atoms sum_orbitals

		println("---------------------------\n")


		local Op  = Operators.Operator(R; sum_atoms=sum_atoms, sum_orbitals=sum_orbitals, kwargs...)


		po(Op)

		v = Op(P) 

		@show size(v)

		@assert size(v)[end]==nr_wf 

		ndims(v)==1 ? push!(X,v) : push!(X,sum(v,dims=1)[:])

	end 

	if length(X)>1  
		
		@assert all(isapprox.(X[1:1],X[2:end]))

		println("\n*************** Test passed\n")

	end 

end 
end 



