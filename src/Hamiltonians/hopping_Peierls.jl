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

function PeierlsPhaseFactor(phi::Real, uc_area::Real, args...
														)::Union{ComplexF64,Function}

	PeierlsPhaseFactor(phi/uc_area, args...)

end


function PeierlsPhaseFactor(B::Real, 
														transl::Union{AbstractString,Char},
														args...
														)::Union{ComplexF64,Function}

	PeierlsPhaseFactor(B, Symbol(lowercase(translation)), args...)

end 



function PeierlsPhaseFactor(B::Real, translation::Symbol,
														)::Union{ComplexF64,Function}

	isapprox(B,0,atol=1e-14) && return 1.0 + 0.0im

	vt = Val(translation)
	
	e = exp(0.5im*pi*B)

	return function p_p_f(ri::AbstractVector,rj::AbstractVector)::ComplexF64

		e^PeierlsIntegral(ri, rj, vt)

	end 

end 



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


Peierls_hopping = HoppingTerm(HamiltBasis(false,false),
															PeierlsPhaseFactor,
															1)




#v0, cond, val, nmax, dmax = Peierls_hopping(params...)

#		Float64
#		Bool/Function 
#		xxxxxxxxx
#		Int 
#		Float64















