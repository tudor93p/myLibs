function swave_gapfunction!(
														ri::AbstractVector{<:Real}, 
													 rj::AbstractVector{<:Real},
													 psi::Number, 
													 args...;
													 kwargs...
													 )::Nothing 

	swave_gapfunction!(ri-rj, psi, args...; kwargs...)

end 


function swave_gapfunction!(
														ri::AbstractVector{<:Real},
													 rj::AbstractVector{<:Real},
													 psi::Function, args...;
													 kwargs...
													 )::Nothing 

	swave_gapfunction!(ri-rj, psi(ri/2+rj/2), args...; kwargs...)

end 


function swave_gapfunction!(
														R::AbstractVector{<:Real},
														psi::Number,
														D::AbstractMatrix{ComplexF64},
														desired_basis::HamiltBasis,
														delta::Function
														)::Nothing  

	all(delta, view(R,1:2)) || return  

	swave_gapfunction!(D, psi, 
										 Val(desired_basis.use_spin), 
										 Val(desired_basis.use_Nambu))

end 


function swave_gapfunction!(D::AbstractMatrix{ComplexF64}, psi::Number,
														 Spin::Val{true}, Nambu::Val{true}
														 )::Nothing

	swave_gapfunction!(view(D, 1:2,3:4), psi, Spin, Val(false)) 

	swave_gapfunction!(view(D, 3:4, 1:2), -conj(psi), Spin, Val(false))
	
end 



function swave_gapfunction!(D::AbstractMatrix{ComplexF64}, psi::Number,
														 Spin::Val{true}, Nambu::Val{false}
														 )::Nothing

	D[1,2] += psi  

	D[2,1] -= psi 

	return 
	
end 

function swave_gapfunction!(D::AbstractMatrix{ComplexF64}, psi::Number,
														 Spin::Val{false}, Nambu::Val{true} 
														 )::Nothing

	D[1,2] += psi
		
	D[2,1] += conj(psi) 

	return 
	
end 

#
#function swave_gapfunction(x::AbstractVector{<:Real},
#													 args...; kwargs...)::Matrix{ComplexF64}
#
#	D = zeros(ComplexF64, 2, 2) 
#
#	swave_gapfunction!(D, x, args...; kwargs...)
#
#	return D
#
#end 
#






swave = HoppingTerm(HamiltBasis(false,true),
										swave_gapfunction!, 0)







