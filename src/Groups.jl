module Groups 
############################################################################# 


import ..LA 




#===========================================================================#
#
# Pauli 2x2 matrices 
#
#---------------------------------------------------------------------------#



function PauliMatrices()::Dict{Int64,Matrix{ComplexF64}}

	Dict(i=> PauliMatrix(i) for i=0:3)

end


PauliMatrix(i::Int)::Matrix{ComplexF64} = PauliMatrix(Val(i))

PauliMatrix(::Val{0})::Matrix{ComplexF64} = [[1 0]; [0 1]]
PauliMatrix(::Val{1})::Matrix{ComplexF64} = [[0 1]; [1 0]]
PauliMatrix(::Val{2})::Matrix{ComplexF64} = [[0 -1im]; [1im 0]]
PauliMatrix(::Val{3})::Matrix{ComplexF64} = [[1 0]; [0 -1]]



#===========================================================================#
#
#	Gamma 4x4 matrices 
#
#---------------------------------------------------------------------------#

function GammaMatrix!(G::AbstractMatrix{ComplexF64},
										 i::Int,j::Int)::Nothing 

	kron!(G,PauliMatrix(i),PauliMatrix(j))

	return 

end 
function GammaMatrix(i::Int,j::Int)::Matrix{ComplexF64}

	kron(PauliMatrix(i),PauliMatrix(j))

end 

#===========================================================================#
#
# general spin Matrices
#
#---------------------------------------------------------------------------#



function SpinMatrices(s::Real)::Dict
	
    """
    Construct spin-s matrices for any half-integer spin.

    Parameters
    ----------

    s : float or int
        Spin representation to use, must be integer or half-integer.
    include_0 : bool (default False)
        If `include_0` is True, S[0] is the identity, indices 1, 2, 3
        correspond to x, y, z. Otherwise indices 0, 1, 2 are x, y, z.

    Returns
    -------

    ndarray
        Sequence of spin-s operators in the standard spin-z basis.
        Array of shape `(3, 2*s + 1, 2*s + 1)`, or if `include_0` is True
        `(4, 2*s + 1, 2*s + 1)`.
    """

	d = Int(2s+1)

	Sz = 1/2 * LA.diagm(0=>d-1:-2:-d)
												
  # first diagonal for general s from en.wikipedia.org/wiki/Spin_(physics)
	
	diag = [1/2*sqrt(i*(d-i)) for i in 1:d-1]

	Sx = LA.diagm(1=>diag,-1=>diag)
	
	Sy = im*LA.diagm(1=>-diag,-1=>diag)

	return Dict(0=>ArrayOps.UnitMatrix(d), 1=>Sx, 2=>Sy, 3=>Sz)

end


#===========================================================================#
#
#  
#
#---------------------------------------------------------------------------#


function WeylVector(M::AbstractMatrix{<:Number})::Vector{<:Number}

	out = [LA.dot(PauliMatrix(j),M)/2 for j=0:3]

	return all(<(1e-14)∘abs∘imag, out) ? real(out) : out 

end 

																												 



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#




function WeylRepr(t::Number, x::Number, y::Number, z::Number)::Matrix{ComplexF64} 

	A  = t*PauliMatrix(0)

	A += x*PauliMatrix(1)

	A += y*PauliMatrix(2)
	
	A += z*PauliMatrix(3)

	return A 

end #generic Hermitian matrix when the args are real 

dWeylRepr(i::Int)::Matrix{ComplexF64} = PauliMatrix(i)  


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


function O2Repr!(R::AbstractMatrix{<:Union{Float64,ComplexF64}}, 
								 theta::Real)::Nothing

	@assert LA.checksquare(R)==2 

	R[1,1] = cos(theta) 

	R[2,2] = R[1,1]

	R[1,2] = sin(theta) 
	
	R[2,1] = -R[1,2]

	return 

end 


function O2Repr(theta::Real)::Matrix{Float64}

	R = zeros(2,2)

	O2Repr!(R, theta)

	return R

end 



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#






#function U2Repr!(out::AbstractMatrix{ComplexF64},
#								 phi::Real, 
#								psi::Real, 
#								theta::Real, 
#								delta::Real 
#								)::Nothing#Matrix{ComplexF64}
#
##	out = zeros(ComplexF64,2,2) 
#
#	@assert LinearAlgebra.checksquare(out)==2 
#
#	O2Repr!(out,theta) 
#
#	out[1,:] *= exp(1im*(phi/2+psi))
#	out[2,:] *= exp(1im*(phi/2-psi)) 
#
#	out[:,1] *= exp(1im*delta)
#	out[:,2] *= exp(-1im*delta)
#
#	return 
#
#end 

function U2Repr(args::AbstractVector{<:Real})::Matrix{ComplexF64}

	U2Repr(args...)

end 

function U2Repr(args::Vararg{<:Real,4})::Matrix{ComplexF64}
	
	out = WeylRepr(args...)

	out .= exp(im*out)

	return out 

end 
















































































































































































































#############################################################################
end 

