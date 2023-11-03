module QuantumMechanics 
#############################################################################

import LinearAlgebra 


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

function sortperm_energy(N::Int;
												 halfspace::Bool=false,
												 occupied::Bool=false,
												 kwargs...
												 )::AbstractVector{Int}

	halfspace || return 1:N 

	return occupied ? (1:div(N,2)) : (div(N,2)+1:N)

end 



function sortperm_energy(E::AbstractVector{<:Real}; 
												 halfspace::Bool=false,
												 kwargs...
												 )::AbstractVector{Int}

	halfspace || return sortperm(E; kwargs...)

	return partialsortperm(E, 
												 sortperm_energy(length(E); 
																				 halfspace=halfspace, kwargs...)
												 )
end 



function psien_sorted_energy(
													 psi::AbstractArray{ComplexF64,N},
														E::AbstractVector{<:Real};
														vsdim::Int=2,
														kwargs...
														)::AbstractVector{<:AbstractArray} where N

	@assert N>=2  

	@assert vsdim<=N 

	i = sortperm_energy(E; kwargs...)

	return [selectdim(psi, vsdim, i), view(E, i)]

end   


function psi_sorted_energy(psi::AbstractArray{ComplexF64,N},
													 E::Union{<:Int,<:AbstractVector{<:Real}};
													vsdim::Int=2,
												 kwargs...
											 )::AbstractArray{ComplexF64,N} where N

	@assert N>=2 

	@assert vsdim<=N

	return selectdim(psi, vsdim, sortperm_energy(E; kwargs...))

end 

function psi_sorted_energy(psi::AbstractArray{ComplexF64,N};
													vsdim::Int=2,
												 kwargs...
											 )::AbstractArray{ComplexF64,N} where N

	psien_sorted_energy(psi, size(psi,vsdim); vsdim=vsdim, kwargs...)

end 


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function WFoverlap(wf1::AbstractMatrix{T1},
								 wf2::AbstractMatrix{T2}=wf1;
								 vsdim::Int=2
								 )::Matrix{promote_type(T1,T2)
													 } where {T1<:Number,T2<:Number} 

	WFoverlap!(
		Matrix{promote_type(T1,T2)}(undef, size(wf1,vsdim), size(wf2,vsdim)),
		wf1, wf2;
		vsdim=vsdim
		)

end 


function WFoverlap!(args...; vsdim::Int=2
								 )::AbstractMatrix{<:Number}

	(vsdim==2 ? WFoverlap2! : WFoverlap1!)(args...)

end 



function WFoverlap2!(A::T,
									wf1::AbstractMatrix{<:Number},
									wf2::AbstractMatrix{<:Number}=wf1,
								 )::T where T<:AbstractMatrix{<:Number}

	LinearAlgebra.mul!(A, wf1', wf2)

end  

function WFoverlap1!(A::T,
									wf1::AbstractMatrix{<:Number},
								 )::T where T<:AbstractMatrix{<:Number}

	WFoverlap2!(A, transpose(wf1))

end 

function WFoverlap1!(A::T,
									wf1::AbstractMatrix{<:Number},
									wf2::AbstractMatrix{<:Number},
								 )::T where T<:AbstractMatrix{<:Number}

	WFoverlap2!(A, transpose(wf1), transpose(wf2))

end 















############################################################################# 
end 

