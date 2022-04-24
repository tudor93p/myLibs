module SignalProcessing 
#############################################################################


import Dierckx, FFTW, ..Algebra 

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


function LinearInterp1D(y1::AbstractArray{T1,N} where T1<:Number, 
												y2::AbstractArray{T2,N} where T2<:Number,
												x1::Real=0, x2::Real=1
												)::Function where N

	x1>x2 && return LinearInterp1D(y2, y1, x2, x1)

	x_to_01 = Interp1D([x1,x2],[0,1],1)

	slope = y2 - y1 

	T = typeof(0.1*sum(slope))

	return function lin_interp_1D(x::Real)::Array{T,N} 

		y1 + slope*x_to_01(x)

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


function identifyPeaksDips_cubicSpline(spline::Dierckx.Spline1D 
																	 )::NTuple{2,Matrix{Float64}}

	@assert spline.k == 3

	knots = Dierckx.get_knots(spline)  
	
	crit_pts = criticalPoints_cubicSpline(spline, knots)   # already sorted 


	D2neg = if isempty(crit_pts) Int[] else 
						
							findall(<(0), Dierckx.derivative(spline, crit_pts, 2))

					end 

	peaks = Matrix{Float64}(undef, 2, length(D2neg))
	dips = Matrix{Float64}(undef, 2, length(crit_pts)-length(D2neg)+2)

  dips[1,1] = knots[1]
  dips[1,2:end-1] = view(crit_pts, setdiff(axes(crit_pts,1), D2neg))
  dips[1,end] = knots[end] 
	dips[2,:] = spline(selectdim(dips,1,1)) 

	
	peaks[1,:] = view(crit_pts, D2neg) 

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

		start = 1 + get(right_bases_int, i-1, 0)

		right_bases_int[i] = findnext(>(x), selectdim(dips,1,1), start)

	end  


	for (i,(n,y)) in enumerate(zip(right_bases_int,selectdim(peaks,1,2)))

		L = findlast(>(y), view(peaks, 2, 1:i-1)) # higher peak left  
		L = (isnothing(L) ? 1 : right_bases_int[L], n-1) 

		R = findfirst(>(y), view(peaks, 2, i+1:size(peaks,2))) # higher right
		R = (n, isnothing(R) ? size(dips,2) : i+right_bases_int[R]-1)

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

























































































































































































































































































#############################################################################
end 


