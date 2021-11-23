

function prep_PeierlsPhaseFactor(phi::Real, uc_area::Real, args...
														)::Tuple

	prep_PeierlsPhaseFactor(phi/uc_area, args...)

end


function prep_PeierlsPhaseFactor(B::Real, 
														transl::Union{AbstractString,Char},
														args...
														)::Tuple

	prep_PeierlsPhaseFactor(B, Symbol(lowercase(translation)), args...)

end 



function prep_PeierlsPhaseFactor(B::Real, translation::Symbol,
														)::Tuple

	isapprox(B,0,atol=1e-14) && return ()

	return (exp(0.5im*pi*B), Val(translation))

end 


function PeierlsIntegral(A::AbstractVector, B::AbstractVector,
												 ::Val{:x})::Float64

	# Landau gauge A = [-y, 0, 0] 
	

		(A[1]-B[1])*(A[2]+B[2])

end 


function PeierlsIntegral(A::AbstractVector, B::AbstractVector,
												 ::Val{:y})::Float64
	
	#	Landau Gauge A = [0, x, 0]
	
	(A[1]+B[1])*(B[2]-A[2])

end
	

#	return function p_p_f(ri::AbstractVector,rj::AbstractVector)::ComplexF64
#
#		e^PeierlsIntegral(ri, rj, vt)
#
#	end 
#
#end 

function PeierlsPhaseFactor(ri::AbstractVector{<:Real}, 
														rj::AbstractVector{<:Real},
														t::Number
														)::ComplexF64

	t	

end 

function PeierlsPhaseFactor(ri::AbstractVector, rj::AbstractVector,
														t::Number,
														e::ComplexF64, vt::Val
													 )::ComplexF64

	t * e^PeierlsIntegral(ri, rj, vt)

end 

#hopp_Peierls.tij(ri, rj, params...)
#isapprox(B,0,atol=1e-14) && return ()
#
#	vt = Val(translation)
#	
#	e = exp(0.5im*pi*B)
#
#	return (e,vt) 
#
#end 





#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


hopping_Peierls = HoppingTerm(HamiltBasis(false,false), PeierlsPhaseFactor, 1)
















