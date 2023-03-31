#using Distributed 
#include("identify_sectors_.jl")


#@show length(Bonds) 


import myLibs: Utils 


#pmap(1:20) do i
#
#	for n in ([100,1000,10000])
#
#Utils.IdentifySectors(first, Bonds[1:n]);
#
#	end 
#		
#	Utils.IdentifySectors(first, Bonds)
#
#
#
#end  

@testset "left/right sectors" begin 


	for t in (Int,Float64), n in 10:15 

		L = rand(t,n)

		A = zeros(t,length(L)+1); 

		cumsum!(view(A,2:length(A)),L);

		@test cumsum([0;L])≈A

	end  

	for n =10:15 

		B = sort(unique(rand(1:n,5)))

		S = Utils.IdentifySectors(B)

		R = last.(S)

		@test Utils.rSectors(first.(S),length(B))==S

		@test Utils.lSectors(R)==S 

		b=cumsum([0;B])

		@test Utils.sepLengths_cumulRanges(B)==[b[i]+1:b[i+1] for i=eachindex(B)]



#		@test Utils.get_sectors(vcat(0,B),1)==Utils.lSectors(B) 
#		@test Utils.get_sectors(vcat(B,length(B)),0)==Utils.rSectors(B) 

	end 

end 


@testset "recursive vs non-recursive" begin 

#rand(1:N,n)

A = [1,1,1,3,4,5,5,1] 


	s1 = Utils.IdentifySectors_customF(==,A)

#	s2 = Utils.IdentifySectors_customF2(==,A)
#
#	@test s1==s2 

	@test Utils.lSectors(last.(s1))==s1

	for n =10:15 

		A = sort(rand(1:n,10))


	s1 = Utils.IdentifySectors_customF(==,A)

#	s2 = Utils.IdentifySectors_customF2(==,A)
#
#	@test s1==s2  

	rB = last.(s1)
	lB = first.(s1)

	@test Utils.lSectors(rB)==s1
	@test Utils.rSectors(lB,length(A))==s1

	q = falses(length(A))

	for i in rB 
		q[i]=1
	end 

	@test s1==Utils.lSectors(q) 

	q.= 0 
	for i in lB
		q[i]=1
	end 
	@test s1==Utils.rSectors(q) 

	end 


	for p =1:0.2:5

		n = Int(round(10.0^p ))

	#	n<50000 || break 

		@show n 

		A = sort!(rand(1:n,n))


@time	s1 = Utils.IdentifySectors_customF(==,A)

#	@time s2 = Utils.IdentifySectors_customF2(==,A)
#@test s1==s2 
println()
end 

end 





