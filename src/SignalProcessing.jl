module SignalProcessing 
#############################################################################


import Dierckx, FFTW, LinearAlgebra 
import ..Algebra

#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


#	dist = Utils.Unique(diff(sort(x)), tol=1e-6)
#
#	length(dist)>1 && return fourier_abs(interp(x, y; dim=dim)...; dim=dim)
#
#	return (
#
#				 FFTW.rfftfreq(length(x), 2pi/(x[2]-x[1])),
#			
#				 abs.(mapslices(FFTW.rfft, y, dims=dim))
#


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#





# 	!!! 'dim' here has the opposite meaning in Lattices !!!
#like mapslices vs eachslice



function fft(x::AbstractMatrix{<:Real}, w::Real=0; 
						 dim::Int, addup=true)::Array{ComplexF64}
	
	Y = ArrayOps.multiply_elementwise(exp.(-1im*w*(axes(x,dim).-1)), x, dim)

	return addup ? dropdims(sum(Y, dims=dim), dims=dim) : Y
	
end 


function fft(x::AbstractVector{<:Real}; kwargs...)::Array{ComplexF64}
	
	fft(x, 2pi*(axes(x,1).-1)/length(x); kwargs...)

end 


function fft(x::AbstractVector{<:Real}, 
						 w::Union{Real,AbstractVector{<:Real}};
						 addup::Bool=true, kwargs...)::Array{ComplexF64}

	E = exp.(-1im*OuterBinary(w, axes(x,1).-1, *))

	return addup ? E*x : E .* reshape(x,1,:)

end 



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

#function Interp1D(x, y, k::Int)::Union{Function,Dierckx.Spline1D}
#
#	if k==0
#
#		get_ind = Interp1D(vcat(x...), 1:length(x), 1)
#
#		out(X::Real) = y[Int(round(get_ind(X)))]
#
#		out(X::AbstractVector) = y[Int.(round.(get_ind(X)))]
#		
#		return out 
#
#	end 
#
#
#	return Dierckx.Spline1D(vcat(x...), vcat(y...), k=k)
#
#end 

function Interp1D_k0(x, y)::Function 

	get_ind = Interp1D(vcat(x...), 1:length(x), 1)

	T = eltype(first(y))

	out(X::Real)::T = y[Int(round(get_ind(X)))]

	out(X::AbstractVector{<:Real})::Vector{T} = [y[Int(round(i))] for i in get_ind(X)]
	
	return out 

end 

function knownerr1(E)::Bool

	isa(E,ErrorException) && E.msg=="The maximal number of iterations maxit (set to 20 by the program)\nallowed for finding a smoothing spline with fp=s has been reached: s\ntoo small. There is an approximation returned but the corresponding\nweighted sum of squared residuals does not satisfy the condition\nabs(fp-s)/s < tol."

end 

function Interp1D_knon0(x::AbstractVector{<:Real}, 
												y::AbstractVector{<:Real}, k::Int;
												s::Real=0.0,
												)::Dierckx.Spline1D 
	try 

		return Dierckx.Spline1D(x, y; k=k, s=s)

	catch E 

		knownerr1(E) || throw(E) 

		for d in 0.01:0.01:1 
			
			for s_ in round.(s .+ [1,-1].*(d .- rand(2)*0.01),digits=4)

				0<= s_ <=1 || continue 
	
				try 
				
					spl = Dierckx.Spline1D(x, y; k=k, s=s_)  
	
					@warn "Interpolation with smoothness $s failed. Used $s_ instead"
	
					return spl 
	
				catch E1 
	
					knownerr1(E1) || throw(E1) 
	
				end 

			end 		
			
		end 

	end 

end 






function Interp1D_knon0(x::AbstractVector{<:Real}, 
												y::AbstractVector{<:ComplexF64}, k::Int;
												s::Real=0.0,
												)::Function 

	f_re = Dierckx.Spline1D(x, real(y); k=k, s=s)
	f_im = Dierckx.Spline1D(x, imag(y); k=k, s=s)

	f(X::AbstractVector{<:Real})::Vector{ComplexF64} = f_re(X) + 1im*f_im(X)

	f(X::Real)::ComplexF64 = f_re(X) + 1im*f_im(X)

	return f 

end 

function Interp1D_knon0(x, y, k::Int;
											 kwargs...)::Union{Function,Dierckx.Spline1D}

	Interp1D_knon0(vcat(x...), vcat(y...), k; kwargs...)

end 

function Interp1D(x, y, k::Int; kwargs...)

	k==0 ? Interp1D_k0(x, y) : Interp1D_knon0(x, y, k; kwargs...)

end 





function Interp1D(x, y, k::Int, X; kwargs...)

	Interp1D(x,y,k; kwargs...)(X)

end


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function LinearInterp1D(y1::T1, y2::T2 
												)::Function where {T1<:Number,T2<:Number}

	function lin_interp_1D(x::Real)::promote_type(T1,T2,Float64)

		x<=0 && return y1 
		x>=1 && return y2 

		return y1 + (y2-y1)*x 

	end 

end 


function LinearInterp1D(y1, y2, x1::Real, x2::Real)::Function 

	x1>x2 && return LinearInterp1D(y2, y1, x2, x1)

	f = LinearInterp1D(y1, y2) 

	x1==0 && x2==1 && return f

	return f ∘ Interp1D([x1,x2],[0,1],1)


end 


function LinearInterp1D(y1::AbstractArray{T1,N},
												y2::AbstractArray{T2,N}
												)::Function where {N,T1<:Number,T2<:Number}

	slope = y2 - y1 


	return function lin_interp_1D(x::Real)::Array{promote_type(Float64,T1,T2),N}

		out = Array{promote_type(Float64,T1,T2),N}(undef, size(slope))

		x<=0 && return copy!(out,y1)
		x>=1 && return copy!(out,y2)

		return LinearAlgebra.axpy!(x, slope, copy!(out, y1))

#		x*slope + out 

#		y1 + slope*x 

	end 

end 


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#



function criticalPoints_cubicSpline(spline::Dierckx.Spline1D,
																		knots::AbstractVector{<:Real} = Dierckx.get_knots(spline)
																	 )::Vector{Float64}

	@assert spline.k==3 
	
	
	
	solutions = Matrix{Float64}(undef, 3, length(knots)-1)
	
	aux_storage = Vector{Float64}(undef, 3) 
	
	three_knots = Vector{Float64}(undef, 3)  
	
	three_values = Vector{Float64}(undef, 3) 
	
	
	three_knots[1] = knots[1]
	
	
	for (ik,k) in enumerate(knots[2:end])
	
		three_knots[3] = k
		
		three_knots[2] = three_knots[1]/2 + three_knots[3]/2
		
		
		for j in UnitRange(1 + (ik>1), 3)
		
			three_values[j] = Dierckx.derivative(spline, three_knots[j]) 
		
		end 
		
		
		Algebra.poly2roots_from3vals!(selectdim(solutions, 2, ik),
													aux_storage, 
													three_knots, 
													three_values, 
													knots[1]-1)
		
		
		three_values[1] = three_values[3]  

		three_knots[1] = three_knots[3] 
	
	end  
	
	return sort(solutions[solutions.>=knots[1]])
	
end 







#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function identifyPeaksDips_cubicSpline(x::AbstractVector{<:Real},
																		 y::AbstractVector{<:Real};
																		 kwargs...
																	 )::NTuple{2,Matrix{Float64}}

	identifyPeaksDips_cubicSpline(Interp1D(x, y, 3; kwargs...))

end  



function findMinMax_criticalPoints(spline::Dierckx.Spline1D,
																	 crit_pts::AbstractVector{<:Real},
																	 )::NTuple{2,BitVector}

	isempty(crit_pts) && return (Int[],Int[])

	derivs = Dierckx.derivative(spline, crit_pts, 2)

	i_max = falses(length(crit_pts))

	i_min = falses(length(crit_pts))

	for (i,d) in enumerate(derivs)

		if d>=0

			i_min[i] = true 

		elseif i==1 || !i_max[i-1]

			i_max[i] = true 

		end 

	end 

	return i_min,i_max 

end 


function getMinMax_criticalPoints(spline::Dierckx.Spline1D,
																	knots = Dierckx.get_knots(spline)  
																	 )::Tuple{<:AbstractVector{<:Real},
																						<:AbstractVector{<:Real}
																						}

	crit_pts = criticalPoints_cubicSpline(spline, knots)   # already sorted  

	i_dips, i_peaks = findMinMax_criticalPoints(spline, crit_pts)

	return view(crit_pts, i_dips), view(crit_pts, i_peaks)

end 


function identifyPeaksDips_cubicSpline(spline::Dierckx.Spline1D 
																	 )::NTuple{2,Matrix{Float64}}

	@assert spline.k == 3

	
	knots = Dierckx.get_knots(spline)  


	x_dips, x_peaks = getMinMax_criticalPoints(spline, knots)


	peaks = Matrix{Float64}(undef, 2, length(x_peaks))
	dips = Matrix{Float64}(undef, 2, length(x_dips)+2)


  dips[1,1] = knots[1]
	dips[1,2:end-1] = x_dips
  dips[1,end] = knots[end] 

	dips[2,:] = spline(selectdim(dips,1,1)) 


	peaks[1,:] = x_peaks

	isempty(peaks) || setindex!(peaks, spline(selectdim(peaks,1,1)), 2, :)

	return peaks, dips 

end 








#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#
function peakProminences((peaks,dips)::Tuple{AbstractMatrix{<:Real},
																						 AbstractMatrix{<:Real}}
												 )::NTuple{3,Matrix{Float64}} 

	peakProminences(peaks, dips)

end 




function peakProminences(peaks::AbstractMatrix{<:Real},
												 dips::AbstractMatrix{<:Real}
												 )::NTuple{3,Matrix{Float64}}

	peaks_proms = vcat(copy(peaks), copy(selectdim(peaks,1,2:2)))

	right_bases_int = zeros(Int, size(peaks,2))

	for (i,x) in enumerate(selectdim(peaks,1,1))

		start = max(1, get(right_bases_int, i-1, 0))

		right_bases_int[i] = findnext(>(x), selectdim(dips,1,1), start)

	end  

	
	for (i,(n,y)) in enumerate(zip(right_bases_int,selectdim(peaks,1,2)))

		L = findlast(>(y), view(peaks, 2, 1:i-1)) # higher peak left  
		L = (isnothing(L) ? 1 : right_bases_int[L]), n-1

		R = findfirst(>(y), view(peaks, 2, i+1:size(peaks,2))) # higher right
#		R = n, (isnothing(R) ? size(dips,2) : i+right_bases_int[R]-1)
		R = n, (isnothing(R) ? size(dips,2) : right_bases_int[R+i-1])

		peaks_proms[3,i] -= maximum(minimum(view(dips,2,a:b)) for (a,b) in (L,R))
		
	end 


	return (peaks_proms, 
					selectdim(dips, 2, right_bases_int.-1),
					selectdim(dips, 2, right_bases_int)
					)

end 







#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


#peakProminences_cubicSpline = peakProminences∘identifyPeaksDips_cubicSpline


function peakProminences_cubicSpline(x::AbstractVector{<:Real},
																		 y::AbstractVector{<:Real};
																		 kwargs...)::NTuple{3,Matrix{Float64}} 

	peakProminences(identifyPeaksDips_cubicSpline(x,y;kwargs...))

end 



function peakProminences_cubicSpline(spl3::Dierckx.Spline1D
																		 )::NTuple{3,Matrix{Float64}} 

	@assert spl3.k==3 

	peakProminences(identifyPeaksDips_cubicSpline(spl3))

end 







#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

# algo: https://www.johndcook.com/blog/standard_deviation/
# 1962 Welford  
# #
# used in: SnakeStates/src/DeviceWithLeads/src/TasksPlots.jl
# and in: Luxemburg/check_WLO/WLO/src/WLO.jl 
# 
#

function init_run_ms()::Vector{Float64}

	zeros(4)

end 

function init_run_m()::Vector{Float64}

	zeros(2)

end 

function init_run_ms!(storage::AbstractVector{Float64})::Nothing 

	storage .= 0 

	return 

end  

function run_m1!(storage::AbstractArray{<:Number}, 
								 k::Real, xk::Number,
								 inds...)::Nothing 

	storage[inds...] += (xk-storage[inds...])/k 

	return 

end 

function run_m1!(storage::AbstractArray{<:Number}, 
								 k::Real, xk::AbstractArray,
								 inds...)::Nothing 

	LinearAlgebra.axpby!(1/k, xk, 1-1/k, view(storage, inds...))

	return 

end 


function run_m!(storage::AbstractVector{Float64}, xk::Real)::Nothing

	storage[1] += 1

	run_m1!(storage, storage[1], xk, 2)

end 


function run_ms!(storage::AbstractVector{Float64}, xk::Real)::Nothing

	# storage[1]: k-1 
	
	storage[1] += 1

	# storage[1]: k

	#	storage[2] and storage[3]: M(k-1)

	storage[3] += (xk-storage[3])/storage[1]

	# storage[3]: M(k) 

	#	storage[4] S(k-1) 
	
	storage[4] += (xk-storage[2])*(xk-storage[3])

	# storage[4]: S(k) 
	
	storage[2] = storage[3] 

	# storage[2]: M(k)
	
	return 

end 

function run_var(storage::AbstractVector{Float64})::Float64 

	storage[1]>1 ? storage[4]/(storage[1]-1) : 0.0

end 


function run_mean(storage::AbstractVector{Float64})::Float64 

	storage[2] 

end 

function run_sigma(storage::AbstractVector{Float64})::Float64 
	
	sqrt(run_var(storage))

end  


function run_ms!(storage::AbstractVector{Float64},
								 X::AbstractVector{<:Real})::Nothing 

	for x in X 

		run_ms!(storage, x) 

	end  

end  


function run_ms(X::AbstractVector{<:Real})::Vector{Float64}

	storage = init_run_ms()

	run_ms!(storage, X)

	return storage

end 
























































































































































































































































































#############################################################################
end 


