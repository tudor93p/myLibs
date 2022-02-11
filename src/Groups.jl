module Groups 
############################################################################# 



import ..Algebra 





#===========================================================================#
#
# 
#
#---------------------------------------------------------------------------#


function WeylRepr(t::Number, x::Number, y::Number, z::Number)::Matrix{ComplexF64} 

	A  = t*Algebra.PauliMatrix(0)

	A += x*Algebra.PauliMatrix(1)

	A += y*Algebra.PauliMatrix(2)
	
	A += z*Algebra.PauliMatrix(3)

	return A 

end 	

dWeylRepr(i::Int)::Matrix{ComplexF64} = Algebra.PauliMatrix(i)  


function dWeylRepr(::Number, ::Number, ::Number,
									 ::Vararg{<:Number,N}
									 )::Vector{Matrix{ComplexF64}} where N

	@assert N<=1 

	dWeylRepr.(1-N:3)

end 	


function WeylRepr(x::Number, y::Number, z::Number)::Matrix{ComplexF64} 

	WeylRepr(0, x, y, z)

end 


function WeylRepr(x::AbstractVector{<:Number})::Matrix{ComplexF64}

	WeylRepr(x...)

end 

function dWeylRepr(x::AbstractVector{<:Number})::Vector{Matrix{ComplexF64}}

	dWeylRepr(x...)

end 



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#



function SU2Repr(theta::Real, n::AbstractVector{<:Real})::Matrix{ComplexF64}

	exp(WeylRepr(0.5im * theta * n))

end 


function dSU2Repr(theta::Real, n::AbstractVector{<:Real}
									)::Vector{Matrix{ComplexF64}}

	x,y,z = [SU2Repr(theta, n)] .* dWeylRepr(0.5im * theta * n)

	return 0.5im*[n[1]*x + n[2]*y + n[3]*z, theta*x, theta*y, theta*z]

end 




function SU2Repr(theta::Real, i::Int)::Matrix{ComplexF64}

	SU2Repr(theta,  setindex!([0,0,0],1,i))

end  


function dSU2Repr(theta::Real, i::Int)::Matrix{ComplexF64}

	0.5im * SU2Repr(theta, i) * dWeylRepr(i)

end  




function SU2Repr(phi::Real, theta::Real, psi::Real)::Matrix{ComplexF64}

	SU2Repr(phi, 3)*SU2Repr(theta, 2)*SU2Repr(psi, 3)

end   

function SU2Repr(ptp::AbstractVector{<:Real})::Matrix{ComplexF64}

	@assert length(ptp)==3

	SU2Repr(ptp...)

end 

function dSU2Repr(phi::Real, theta::Real, psi::Real
								 )::Vector{Matrix{ComplexF64}}
	
	a = SU2Repr(phi, 3)
	b = SU2Repr(theta, 2) 
	c = SU2Repr(psi, 3)

	return 0.5im*[a*dWeylRepr(3)*b*c, 
								a*b*dWeylRepr(2)*c, 
								a*b*c*dWeylRepr(3)]


end  

function dSU2Repr(ptp::AbstractVector{<:Real})::Vector{Matrix{ComplexF64}}

	@assert length(ptp)==3

	dSU2Repr(ptp...)

end 





#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#




























































































































































































































#############################################################################
end 

