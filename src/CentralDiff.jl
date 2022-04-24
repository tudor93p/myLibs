module CentralDiff 
#############################################################################

import ..Utils, ..LA 




#===========================================================================#
#
# coefficients and weights
#
#---------------------------------------------------------------------------#

const CD_P = (CartesianIndices((0:1,)), CartesianIndices((0:1,0:1)))
const CD_S = (hcat([1,1], [-1,1]),
							hcat([1,1,1,1], [-1,1,-1,1,], [-1,-1,1,1,]))

function central_diff_PS(h::Real,s::Real
											)::Tuple{Matrix{Int},Matrix{Int}}

	(
	hcat([0,0],	 [1,0], [0,1], [1,1]), 
	 hcat([1,1,1,1], [-1,1,-1,1,], [-1,-1,1,1,])
	 )

end  

function central_diff_PS(h::Real
											)::Tuple{Matrix{Int},Matrix{Int}}

	(hcat([0], [1]), hcat([1,1], [-1,1]))

end  

function central_diff_w(args::Vararg{Real,N})::Matrix{Float64} where N

	#	[0.25,0.5/h,0.5/s] 
	
	hcat(0.5/N, (1/(N*h) for h in args)...)  # division by N? 

end 

function central_diff_PW1(args...
												 )::Tuple{Matrix{Int},Matrix{Float64}} 
	
	P,S = central_diff_PS(args...)

	return P, central_diff_w(args...).*S 

end 

function central_diff_PW2(args...)::Tuple{Matrix{Int},Matrix{Float64}}
	
	P,S = central_diff_PS(args...)

	return P, volume_element(args...)*S 

end 

function volume_element(steps::Vararg{Real})::Float64 

	prod(steps; init=1.0)

end 


#===========================================================================#
#
# instructions to use P 
#
#---------------------------------------------------------------------------#

function xyz_neighb( I::CartesianIndex,
											 P::AbstractMatrix{Int},
											 k::Int)#::Tuple{Vararg{Int}}
			
	CartesianIndex(Tuple(i+P[n,k] for (n,i) in enumerate(Tuple(I))))
				
end 


function xyz_neighb( I::NTuple{N,Int},
											 P::AbstractMatrix{Int},
											 k::Int)::NTuple{N,Int} where N#{T<:Number,N}
										
	Tuple(I[n]+P[n,k] for n=1:N) 

end 



#===========================================================================#
#
# iter utils 
#
#---------------------------------------------------------------------------#




#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function collect_midval_deriv_1D(
													A::AbstractArray{T,N1}, 
													steps::Vararg{Real,N}
													)::Array{promote_type(T,Float64),N
																	 } where {T<:Number,N,N1}

	@assert N in 1:2 
	@assert N1==N+1 

	D = Array{promote_type(T,Float64), N}(undef, (size(A)[2:N1].+1)...) 

	collect_midval_deriv_1D!(D, A, steps...)

	return D 

end 

function collect_midval_deriv_1D!(
													D::AbstractArray{T,N},
													A::AbstractArray{<:Number,N1}, 
													steps::Vararg{Real,N}
													)::Nothing where {N,N1,T<:Union{Float64,ComplexF64}}

	@assert N in 1:2 
	@assert N1==N+1 
	@assert size(D) == size(A)[2:N1] .+ 1

	D .= T(0.0)

	dv = volume_element(steps...) 

	@simd for xyz=CartesianIndices(axes(A)[2:N1])
		@simd	for a=1:N+1
			@simd for ngh=1:2^N
				
				@inbounds D[xyz + CD_P[N][ngh]] += CD_S[N][ngh, a] * A[a, xyz] * dv

			end 
		end 
	end 

	return 

end 

function collect_midval_deriv_2D!(
								D::AbstractMatrix{T},
								A::AbstractArray{<:Number,N2},
								 steps::Vararg{Real,N}
								 )::Nothing  where {T<:Union{Float64,ComplexF64},N,N2}

	@assert N in 1:2 
	@assert N2==N+2 
	
	s = size(A)[3:N2]  

	Li = LinearIndices(s.+1) # trivial when N=1 

	@assert LA.checksquare(D)==length(Li) 

	D .= T(0)

	dv = volume_element(steps...)

	@simd for xyz=CartesianIndices(s)
		@simd for a2=1:N+1
			@simd for a1=1:N+1
				@simd for ngh2=1:2^N
					@simd for ngh1=1:2^N

						@inbounds D[Li[xyz + CD_P[N][ngh1]], 
												Li[xyz + CD_P[N][ngh2]] 
												] += *(CD_S[N][ngh1, a1],
															 CD_S[N][ngh2, a2],
															 A[a1,a2, xyz],
															 dv
															 )
					end 
				end 
			end 
		end 
	end 

	return 

end 

function collect_midval_deriv_2D(A::AbstractArray{T,N2},
								 steps::Vararg{Real,N}
								 )::Matrix{promote_type(T,Float64)} where {T<:Number,N,N2}
	
	@assert N in 1:2 
	@assert N2==N+2 

	
	n = prod(size(A)[3:N2].+1)
	
	D = Matrix{promote_type(T,Float64)}(undef, n, n)

	collect_midval_deriv_2D!(D, A, steps...)

	return D

end 


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#





function midval_and_deriv(g::AbstractArray{T,N},
													steps::Vararg{Real,N}
													)::Array{promote_type(T,Float64),N+1
																						 } where {T<:Number,N}
	@assert N in 1:2 
	
	A = Array{promote_type(T,Float64), N+1}(undef, N+1, (size(g).-1)...)

	midval_and_deriv!(A, g, steps...)

	return A 

end 

function midval_and_deriv!(A::AbstractArray{T,N1},
													 g::AbstractArray{T,N},
													steps::Vararg{Real,N}
													)::Nothing where {T<:Number,N,N1}

	@assert N in 1:2 
	@assert N1==N+1 

	s = size(g) .- 1

	@assert size(A)==(N+1, s...)
	
	A .= T(0) 

	w = central_diff_w(steps...) 

	@simd for xyz=CartesianIndices(s) 
		@simd for a=1:N+1 
			@simd for ngh=1:2^N 

				@inbounds A[a, xyz] += g[xyz + CD_P[N][ngh]] * CD_S[N][ngh,a] * w[a]

			end 
		end 
	end 
	
	return 

end 








#===========================================================================#
#
# Two types of containers for the middle values and differences 
#
#---------------------------------------------------------------------------#


function mvd_container(MXY::AbstractArray)::Base.Generator 
	eachslice(MXY, dims=1)

end  

function mvd_container(MXY::AbstractArray{T},
																	 I::Tuple{Vararg{Int}}
																	 )::AbstractVector{T}  where T
	view(MXY, :, I...)

end  

function mvd_container(MXY::AbstractArray{T},
																	 I::CartesianIndex
																	 )::AbstractVector{T}  where T
	view(MXY, :, I)

end  

function mvd_container(MXY::T)::T where T<:Tuple{Vararg{AbstractArray}}

	MXY 

end 

function mvd_container(MXY::T, I::Tuple{Vararg{Int}}
																	 )::Base.Generator where T<:Tuple{Vararg{AbstractArray}}

	(M[I...] for M in MXY)

end 


function mvd_container(MXY::T, I::CartesianIndex 
																	 )::Base.Generator where T<:Tuple{Vararg{AbstractArray}}

	(M[I] for M in MXY)

end 






#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function mid_Riemann_sum(data,
												 F::Function,
												 A::Union{AbstractArray{T,N1},
																	Tuple{N1,AbstractArray{T,N}}},
												 steps::Vararg{Real,N}
												 )::promote_type(T,Float64) where {T<:Number,N,N1}

	@assert N+1==N1 
	@assert N in 1:2 

	out::promote_type(T,Float64) = 0.0 

	for I in CartesianIndices(first(mvd_container(A)))
	
		out += F(data, mvd_container(A, I)...)

	end 

	return out*volume_element(steps...)

end  


function mid_Riemann_sum(data,
												 F::Function,
												 A::Union{AbstractArray{T,N1},
																	Tuple{N1,AbstractArray{T,N}}},
												 output_size::NTuple{M,Int},
												 steps::Vararg{Real,N}
												 )::Array{promote_type(T,Float64),M
																	} where {T<:Number,N,M,N1}

	@assert N+1==N1 
	@assert N in 1:2 
	@assert M>0

	out = zeros(promote_type(T,Float64), output_size...)

	for I in CartesianIndices(first(mvd_container(A)))
	
		out += F(data, mvd_container(A, I)...)

	end 

	return out*volume_element(steps...)

end  





#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#



function eval_fct_on_mvd(data, 
												F::Function,
												A::Union{AbstractArray{T,N1},
																		NTuple{N1,AbstractArray{T,N}}},
												args...
												)::Array where {T<:Number,N,N1}

	@assert N+1==N1 
	@assert N in 1:2 
	
	
	out = init_array_fct_on_mvd(mvd, output_size, args...)

	s = size(first(mvd_container(A)))

	I0 = fill(Colon(),ndims(out)-length(s))

	for I in CartesianIndices(s)

		out[I0..., I] .= F(data, mvd_container(A,I)...)

	end 
	
	return out 

end 


function eval_fct_on_mvd!(
												out::AbstractArray{T1,MN},
												data, 
												F!::Function,
												A::Union{AbstractArray{T2,N1}, 
																 NTuple{N1,AbstractArray{T2,N}}},
												output_size::NTuple{M,Int},
												steps::Vararg{Real,N}
												)::Nothing  where {T1<:Number,T2<:Number,N,N1,M,MN}

	@assert N+1==N1 
	@assert N in 1:2 
	@assert MN==M+N 

	out .= T1(0)

	I0 = fill(Colon(),M)

	for I in CartesianIndices(axes(out)[M+1:MN])

		F!(view(out, I0..., I), data, mvd_container(A,I)...)

	end 
	
	return 

end 









#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function eval_deriv_on_mvd!(A::AbstractArray{T,N12} where T<:Number,
														 data, F!::Function,
														 mvd::AbstractArray{<:Real,N1}, 
														 steps::Vararg{Real,N}
														 )::Nothing where {N,N1,N12}


	@assert N in 1:2 
	@assert N1==N+1 
	
	order = N12-N 
	@assert order in 1:2 

	eval_fct_on_mvd!(A, data, F!, mvd, ntuple(i->N1,order), steps...)

	for (j,w) in enumerate(central_diff_w(steps...)), d=1:order 
	
		selectdim(A,d,j) .*= w

	end 

	return  

end  



function init_array_fct_on_mvd( 
												desired_type::DataType,
												output_size::NTuple{M,Int},
												mvd::Union{AbstractArray{T,N1},
																	 NTuple{N1,AbstractArray{T,N}}},
												steps::Vararg{Real,N}
												)::Array{promote_type(desired_type,T,Float64),N+M
																 } where {T<:Number,N,N1,M}

	@assert N+1==N1 
	@assert N in 1:2  

	Array{promote_type(desired_type,T,Float64), N+M
				}(undef,
					output_size..., 
					size(first(mvd_container(mvd)))...)

end 

function init_array_fct_on_mvd(output_size::Union{Int, Tuple{Vararg{Int}}}, 
															 args...
															 )::Array 

	init_array_fct_on_mvd(Float64, output_size, args...)

end 


function init_array_fct_on_mvd( 
												desired_type::DataType,
												output_size::Int,
												mvd::Union{AbstractArray{<:Number,N1},
																	 NTuple{N1,AbstractArray{<:Number,N}}},
												steps::Vararg{Real,N}
												)::Array where {N1,N}

	init_array_fct_on_mvd(desired_type, ntuple(i->N1,output_size), mvd, steps...)

end 






#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function derivs12_on_mvd!(aux::AbstractArray{<:Number,N2},
													dF::AbstractArray{<:Number,N},
													d2F::AbstractMatrix{<:Number},
													Data, df!::Function, d2f!::Function,
													mvd::AbstractArray{<:Real,N1}, 
													steps::Vararg{Real,N}
													)::Nothing where {N,N1,N2}

	eval_deriv_on_mvd!(selectdim(aux, 1, 1), Data, df!, mvd, steps...)

	collect_midval_deriv_1D!(dF, selectdim(aux,1,1), steps...)   

	eval_deriv_on_mvd!(aux, Data, d2f!, mvd, steps...)

	collect_midval_deriv_2D!(d2F, aux, steps...)


end 


function derivs12_on_mvd!(aux::AbstractArray{<:Number,N2},
													Data, df!::Function, d2f!::Function,
													mvd::AbstractArray{<:Real,N1}, 
													steps::Vararg{Real,N}
													)::Tuple{Array{<:Number,N},
																	 Matrix{<:Number}
																	 } where {N,N1,N2}

	eval_deriv_on_mvd!(selectdim(aux, 1, 1), Data, df!, mvd, steps...)

	dF = collect_midval_deriv_1D(selectdim(aux,1,1), steps...)   

	eval_deriv_on_mvd!(aux, Data, d2f!, mvd, steps...)

	d2F = collect_midval_deriv_2D(aux, steps...)

	return (dF,d2F)

end 
	

function derivs12_on_mvd!(Data, df!::Function, d2f!::Function,
													mvd::AbstractArray{<:Real,N1},
													T::DataType,
													steps::Vararg{Real,N}
													)::Tuple{Array{<:Number,N},
																	 Matrix{<:Number},
																	 Array{<:Number,N+2},
																	 } where {N,N1}

	aux = init_array_fct_on_mvd(T, 2, mvd, steps...)

	return (derivs12_on_mvd!(aux, Data, df!, d2f!, mvd, steps...)..., aux)

end 

	









#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function numerical_derivative(f, x, args...
															)::Union{Number,AbstractArray{<:Number}} 

	@assert applicable(f, x)

	return numerical_derivative_(f, x, args...)

end 

function numerical_derivative_(f,
															x::Union{Number,AbstractArray{<:Number}},
															dx::Float64,
															i1::Int, inds::Vararg{Int}
															)::Union{Number,AbstractArray{<:Number}}

	fstep(s) = setindex!(zero(x),s,i1, inds...) 

	return numerical_derivative_(f,x,dx,fstep)

end 





function numerical_derivative_(f,
															x::Union{Number,AbstractArray{<:Number}},
															dx::Float64,
															fstep::Function=identity
															)::Union{Number,AbstractArray{<:Number}}

#	coef2 = [-1/2, 0, 1/2]

#	coef8 = [1/280, -4/105, 1/5, -4/5, 0, 4/5, -1/5, 4/105, -1/280]

	coef = [1/12, -2/3, 2/3, -1/12]

	N = div(length(coef),2) 

	ns = vcat(-N:-1,1:N)


	return sum(a * f(x + fstep(n*dx)) for (n,a) in zip(ns,coef))/dx 

end 



function test_derivative(f, f_truth, x, fstep...
												)::Bool




#	@show length(ns) length(coef )


	truth = f_truth(x) 

#	@show truth 


	orders = map(Utils.logspace(1e-2,1e-9,20)) do dx 

		a = numerical_derivative(f, x, dx, fstep...)

		return -log10.([dx, LA.norm(truth - a)]) 

	end  



	for ords in orders 

		ord0,ord = ords 

		if ord < ord0/2 #&& ord < 3
			
			for item in orders 
		
				println(join(round.(item, digits=1),"\t"))  
		
			end  

			return false 

		end 

	end 

	return true 

end 





































#############################################################################
end 


