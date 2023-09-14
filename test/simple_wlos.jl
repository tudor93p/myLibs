import LinearAlgebra 
import Combinatorics 


import myLibs: Groups,MeshInPlace  

PM = Groups.PauliMatrices()


#function rest_permutations(i1::AbstractVector{<:Int}, N::Int,
#													 )::Combinatorics.Permutations{Vector{Int}}
#
#	n1 = length(i1)
#
#	@assert n1<=n
#
#	if n1==n || !allunique(i1)
#		
#		return Combinatorics.permutations(Int[])
#
#	end 
#
#
#	return Combinatorics.permutations(setdiff(1:N,first_inds))
#
#end 
#
#function single_missing_permut_elem(i1::AbstractVector{<:Int}
#															)::Int 
#
#	@assert allunique(i1)
#	
#	return only(setdiff(1:length(i1)+1,i1))
#
#end 
#
#
#



#function 
#
##%	parity = falses(length(rest_inds))
#
#	I = Matrix{Int}(undef,N,length(I2))
#
#	for (c,i2) in enumerate(I2)  
#
#		setindex!(I, i1, 1:n1, c)
#		setindex!(I, i2, n1+1:N, c) 
#
#	
#	end 
#
#	return I
#
#end 
#
#parity[c] = Combinatorics.levicivita(selectdim(I,2,c))
#A += Combinatorics.levicivita([i,j,k])* coll[k] 































