import myLibs:Operators ,Utils,TBmodel
import LinearAlgebra 
using BenchmarkTools


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

@testset "WF" begin 

for psi in eachslice(P,dims=dim)

	@test isapprox(1,LinearAlgebra.norm(psi))

end  

end 
kwargs = (dim=dim,nr_at=nr_at,nr_orb=nr_orb,size_H=size_H)
#
#
#for Op in [#Operators.Position_Expectation(1, atoms; kwargs...),
#					 #Operators.Operator([-1,0,1]; kwargs...),
#					 #Operators.IPR(;kwargs...),
#					 #Operators.LDOS([1]; kwargs...),
#					 ]
#
#
#	@show Op
#
#
#end 
#
#
#for n in 0:3
#
#	@show n
#
#	ks = Utils.Random_Items([:nr_orb,:nr_at,:size_H],n)
#
#	@show ks 
#
#
#	kw = Utils.dict_diff(kwargs,ks)
#
#	@show kw 
#
#	@show Operators.get_nratorbH(;kw...) 
#
##	if n>0
##
##		k = rand(ks) 
##	
##		@show k 
##		@show Operators.get_nratorbH(kwargs[k];kw...) 
##	
##	end 
#
#
#i = rand(2:length(kwargs))
#
#	@show i 
#
#	k = keys(kwargs)[i]
#
#	@show k
#	v = kwargs[k]
#@show v 
#
#	A = rand(v)
#
#
#					 
#
#	@show size(A)
#
#	@show Operators.Operator(A; kwargs...)
#
#	println()
#	println()
#end  
#
#
#
#
#for k in [:nr_orb,:nr_at,:size_H]
#
#	nrs,enough_data = Operators.get_nratorbH(;kwargs...)
#
#	@show k Operators.which_diag(rand(kwargs[k]),nrs)
#
#end 

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


@testset "Expect, vals." begin 

for n in [1,nr_at,nr_orb,size_H] 
	
for N in [(n,), (n,n)]

	R = round.(rand(N...), digits=2) 

	X = [] 

	for (sum_atoms,sum_orbitals) in [
																	 (true,true);
								(length(N)==1 || n==nr_at) ? [(true,false)] : [];
								(length(N)==1 || n==nr_orb) ? [(false,true)] : [];
							 length(N)==1 ? [(false,false)] : [];
																	 ]
#		println("\n---------------------------")
#
#		@show N sum_atoms sum_orbitals
#
#		println("---------------------------\n")


		local Op  = Operators.Operator(R; sum_atoms=sum_atoms, sum_orbitals=sum_orbitals, kwargs...)


#		po(Op)

		v = Op(P) 

		@test size(v,2)==nr_wf 

		sum_atoms && sum_orbitals && @test size(v,1)==1
		sum_atoms && !sum_orbitals && @test size(v,1)==nr_orb
		!sum_atoms && sum_orbitals && @test size(v,1)==nr_at 
		!sum_atoms && !sum_orbitals && @test size(v,1)==size_H
		

		push!(X,sum(v,dims=1)[:])

	end 

	length(X)>1 && @test all(isapprox.(X[1:1],X[2:end]))


end 
#println()
end 

println()
end  






@testset "Trace" begin


for n in [1,nr_at,nr_orb,size_H], 
	N in [(n,), (n,n)],
	M in [rand(size_H), rand(size_H,size_H)]

		R = round.(rand(N...), digits=2) 

	X =[] 

	for what in ["orbitals","atoms"], sum_up in [true,false]

#		println("\n---------------------------")
#
#		@show N what sum_up
#
#		println("---------------------------\n")

		
		sum_atoms,sum_orbitals = what=="orbitals" ? (sum_up,true) : (true,sum_up)


		tr = Operators.Trace(what, R; sum_up=sum_up, kwargs...)

		v1 = tr(M) 
#		@time tr(M) 

#		tr2 = Operators.Trace2(what, R; sum_up=sum_up, kwargs...)




#		po(Op)

#		v2 = tr2(M)
#		@time tr2(M)

#		@show size(v1) size(v2)

#		@test size(v1)==size(v2)

		for v in [v1]
	
			if v isa AbstractVector 
				push!(X,sum(v))

			elseif v isa Number 
				push!(X,v)
			else

				error()
			end 
		end 

	end 



	length(X)>1 && @test all(isapprox.(X[1:1],X[2:end])) 


end 




println()

end 




@testset "Usual operators" begin 

	atom =rand(1:nr_at)

	wf = rand(1:nr_wf)

	p = P[TBmodel.Hamilt_indices(1:nr_orb, atom, nr_orb),wf]


	ldos = Operators.LDOS(;kwargs...)(P) 
	
	@test size(ldos)==(nr_at,nr_wf)

	@test all(isapprox.(imag(ldos),0,atol=1e-12))


	@test ldos[atom,wf]==sum(abs2,p)
	


	dir = rand(1:2)

	x = Operators.Position(dir, atoms; kwargs...)(P)

	


	@test size(x)==(1,nr_wf,)

	@test all(isapprox.(imag(x),0,atol=1e-12))

	@test isapprox(x[wf],LinearAlgebra.dot(selectdim(atoms,dim2,dir),ldos[:,wf]))





	charge = rand(nr_orb).-0.5 

	c = Operators.SingleAtomCharge(charge, atom; kwargs...)(P)

	@test size(c)==(1,nr_wf,)

	@test isapprox(c[wf],LinearAlgebra.dot(charge,abs2.(p)))


	ipr = Operators.IPR(;kwargs...)(P)


	@test size(ipr)==(1,nr_wf,)

	t = Operators.Trace(:orbitals; kwargs...)(abs2.(P[:,wf]))

	@test isapprox(t,ldos[:,wf])

	@test isapprox(ipr[wf],1/sum(abs2,t))
															 




#	rand(1:nr_at)










end 



