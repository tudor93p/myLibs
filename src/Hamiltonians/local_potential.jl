function total_local_pot(ri::AbstractVector{<:Real},
												 rj::AbstractVector{<:Real},
												 ChemicalPotential::Real,
												 LocalPotential::Real,
												 )::Float64 

	LocalPotential + ChemicalPotential 

end 

function total_local_pot(ri::AbstractVector{<:Real},
												 rj::AbstractVector{<:Real},
												 ChemicalPotential::Real,
												 LocalPotential::Function,
												 )::Float64 

	total_local_pot(ri, rj, 
									ChemicalPotential,
									LocalPotential(ri/2+rj/2), 
									)

end 




#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


local_potential = HoppingTerm(HamiltBasis(false,false), total_local_pot, 0)



















