function args_local_potential(LocalPotential::Number, 
						 ChemicalPotential::Number, 
						 )::Tuple{Float64, Float64}

	(LocalPotential + ChemicalPotential, 1.0)

end 



function args_local_potential(LocalPotential::Function, 
															ChemicalPotential::Number,
															)::Tuple{Float64, Function}

#	mu = TBmodel.matrixval(ChemicalPotential)


	(1.0, 

	 function total_local_pot(ri::AbstractVector, rj::AbstractVector
																	 )::Float64#AbstractMatrix
		 
		 LocalPotential(ri,rj) + ChemicalPotential
	

	 end 

	 )

end 

function args_local_potential(;LocalPotential::Union{Number,Function},
															 ChemicalPotential::Number,
															 kwargs...)

	args_local_potential(LocalPotential, ChemicalPotential)

end 







#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


local_potential = HoppingTerm(HamiltBasis(false,false),
															args_local_potential, 
															0)

























