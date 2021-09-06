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








