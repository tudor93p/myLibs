module QuantumMechanics 
#############################################################################


import ..Algebra 



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function FermiDirac(mu::Real=0, T::Real=0)::Union{Algebra.myDistrib,Function}

	@assert T>=0 

	T<1e-20 && return Algebra.Heaviside(-mu)∘-  

	fd(E::Real)::Float64 = Algebra.FermiDirac(mu, E, T)
	
	fd(E::AbstractVector{<:Real})::Vector{Float64} = Algebra.FermiDirac(mu, E, T)

	return fd 
	
end  


hole_occup(el_occup::Real)::Float64 = 1 - el_occup

hole_occup!(el_occup::Real)::Float64 = 1 - el_occup


function hole_occup(el_occup::AbstractVector{<:Real})::Vector{Float64} 
	
	hole_occup.(el_occup)

end 


function hole_occup!(el_occup::T)::T where T<:AbstractVector{Float64}

	map!(hole_occup, el_occup, el_occup)

end 




function FermiDiracHole(args::Vararg{Real,N})::Function where N 
	
	@assert N<=2 

	hole_occup!∘FermiDirac(args...)

end  



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#
















############################################################################# 
end 

