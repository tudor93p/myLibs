#function args_local_potential(LocalPotential::Number, 
#						 ChemicalPotential::Number, 
#						 )::Tuple{Float64, Function}
#
#	(LocalPotential + ChemicalPotential, TBmodel.constHoppTerm(1))
#
#end 
#
#
#
#function args_local_potential(LocalPotential::Function, 
#															ChemicalPotential::Number,
#															)::Tuple{Float64, Function}
#
##	mu = TBmodel.matrixval(ChemicalPotential)
#
#
#	(1.0, 
#
#	 function total_local_pot(ri::AbstractVector, rj::AbstractVector
#																	 )::Float64#AbstractMatrix
#		 
#		 LocalPotential(ri,rj) + ChemicalPotential
#	
#
#	 end 
#
#	 )
#
#end 
#


function total_local_pot(ri::AbstractVector{<:Real},
												 rj::AbstractVector{<:Real},
												 LocalPotential::Real,
												 ChemicalPotential::Real)::Float64 

	LocalPotential + ChemicalPotential 

end 

function total_local_pot(ri::AbstractVector{<:Real},
												 rj::AbstractVector{<:Real},
												 LocalPotential::Function,
												 ChemicalPotential::Real)::Float64 

	total_local_pot(ri, rj, LocalPotential(ri/2+rj/2), ChemicalPotential)

end 

#function args_local_potential(;LocalPotential::Union{Number,Function},
#															 ChemicalPotential::Number,
#															 kwargs...)
#
#	args_local_potential(LocalPotential, ChemicalPotential)
#
#end 







#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


local_potential = HoppingTerm(HamiltBasis(false,false), total_local_pot, 0)



















