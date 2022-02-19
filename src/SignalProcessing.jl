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


function Interp1D_knon0(x::AbstractVector{<:Real}, 
												y::AbstractVector{<:Real}, k::Int;
												s=0.0,
												)::Dierckx.Spline1D 

	Dierckx.Spline1D(x, y; k=k, s=s)

end 


function Interp1D_knon0(x::AbstractVector{<:Real}, 
												y::AbstractVector{<:ComplexF64}, k::Int;
												s=0.0,
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


function LinearInterp1D(y1::AbstractArray{T1,N}, 
												y2::AbstractArray{T2,N},
												x1::Real=0, x2::Real=1
												)::Function where {T1<:Number,T2<:Number,N} 

	x1>x2 && return LinearInterp1D(y2, y1, x2, x1)

	x_to_01 = Interp1D([x1,x2],[0,1],1)

	slope = y2 - y1 

	T = typeof(0.1*slope[1])

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




function identifyPeaksDips_cubicSpline(spline::Dierckx.Spline1D
																	 )::NTuple{2,Matrix{Float64}}

	@assert spline.k == 3

	knots = Dierckx.get_knots(spline)  
	
	crit_pts = criticalPoints_cubicSpline(spline, knots)   # already sorted 


	D2neg = findall(<(0), Dierckx.derivative(spline, crit_pts, 2))


	peaks = Matrix{Float64}(undef, 2, length(D2neg))
	dips = Matrix{Float64}(undef, 2, length(crit_pts)-length(D2neg)+2)


	peaks[1,:] = view(crit_pts, D2neg) 

	dips[1,1] = knots[1] 
	dips[1,2:end-1] = view(crit_pts, setdiff(axes(crit_pts,1), D2neg))
	dips[1,end] = knots[end]

	peaks[2,:] = spline(selectdim(peaks,1,1))
	dips[2,:] = spline(selectdim(dips,1,1))

	return peaks, dips 

end 








#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#
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
































































































































































































































































































#############################################################################
end 


