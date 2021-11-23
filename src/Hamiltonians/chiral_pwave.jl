

function chiral_pwave_gapfunction!(ri::AbstractVector{<:Real},
																	rj::AbstractVector{<:Real},
																	eta::Function,
																	args...;
																	kwargs...
																	)::Nothing
																	
	chiral_pwave_gapfunction!(ri-rj, eta(ri/2+rj/2), args...; kwargs...)

end


function chiral_pwave_gapfunction!(
																	 ri::AbstractVector{<:Real},
																	rj::AbstractVector{<:Real},
																	eta::AbstractVector{<:Number},
																	args...;
																	kwargs...
																	)::Nothing

	chiral_pwave_gapfunction!(ri-rj, eta, args...; kwargs...)

end 




function chiral_pwave_gapfunction!(
																	 R::AbstractVector{<:Real},
																	 eta::AbstractVector{<:Number},
																	 D::AbstractMatrix{ComplexF64},
																	 desired_basis::HamiltBasis,
																	delta::Function,
																	)::Nothing

	Spin = desired_basis.use_spin ? 2 : 1 


	for ((ri,rj),e) in zip((R[1:2],R[2:-1:1]),eta)

		s = e*(delta(ri-1) - delta(ri+1))*delta(rj)*0.5im

		for (i,j) in zip(1:Spin, 2*Spin:-1:Spin+1)

			D[i,j] += s

			if desired_basis.use_Nambu  

				D[j,i] -= conj(s) 

			end 

		end 

	end 

	return 

	# PyCall.pyimport("get_SCgaps").SCmatr_realspace(dist_tol,hopp_cutoff,funs="meta")
end





#function chiral_pwave_gapfunction(desired_basis::HamiltBasis,
#																	args...; kwargs...
#																	)::Matrix{ComplexF64}
#
#	D = zeros(ComplexF64, desired_basis.matrix_dim, desired_basis.matrix_dim)
#
#	chiral_pwave_gapfunction!(args..., D, desired_basis; kwargs...)
#
#	return D 
#
#end 



chiral_pwave = HoppingTerm(HamiltBasis(false,true), 
													 chiral_pwave_gapfunction!, 1) 




