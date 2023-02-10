module Geometry
#############################################################################

import ConcaveHull 

import LinearAlgebra; const LA=LinearAlgebra

#import GeometryBasics 
import LazySets 


import ..Utils, ..ArrayOps, ..Algebra





#===========================================================================#
#
# Decides whether a point is in a polygon
#
#---------------------------------------------------------------------------#


function is_left(P0, P1, P2)::Real

	(P1[1] - P0[1]) * (P2[2] - P0[2]) - (P2[1] - P0[1]) * (P1[2] - P0[2])

end 

function PointInPolygon_wn(P::Utils.List, V::AbstractMatrix; 
													 dim, prepare_vertices=true, kwargs...)::Bool


	prepare_vertices && return PointInPolygon_wn(P, prepare_polygon_vertices(V; dim=dim, kwargs...); dim=dim, prepare_vertices=false, kwargs...)

  wn = 0   # the winding number counter

  						# repeat the first vertex at the end


	v(i) = selectdim(V, dim, i)

	v(i,j) = v(i)[j]



  # loop through all edges of the polygon

	for i in 1:size(V, dim)-1  # edge from v(i) to v(i+1) 

		if v(i,2) <= P[2] # start y <= P[1]

			if v(i+1,2)>P[2] && is_left(v(i), v(i+1),P)>0 
									# an upward crossing & P left of edge
			
				wn += 1           # have a valid up intersect

			end 

			 
		else # start y > P[2] (no test needed)
			 
			if v(i+1,2)<=P[2] && is_left(v(i), v(i+1), P)<0

						
						# a downward crossing &  P right of edge

				wn -= 1           # have a valid down intersect


			end

		end 

	end 

  return wn != 0


end 




#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function Order_PolygonVertices(V::AbstractMatrix{T}; dim)::Matrix{T} where T<:Number

	C = sum(V, dims=dim)[:]/size(V,dim) # center of the polygon

	return sortslices(V, dims=dim, by = R->atan(R[2]-C[2], R[1]-C[1]))

end 



function prepare_polygon_vertices(V::AbstractMatrix; 
																	order_vertices=false, dim)::Matrix

	if order_vertices

		return prepare_polygon_vertices(
								Order_PolygonVertices(V, dim=dim),
								order_vertices=false, dim=dim)

	end 

	return cat(selectdim(V, [2,1][dim], 1:2),
						 selectdim(selectdim(V, [2,1][dim], 1:2), dim, 1:1),
						 dims=dim)
end 

function PointInPolygon_wn(V::AbstractMatrix; kwargs...)::Function

	V = prepare_polygon_vertices(V; kwargs...) 
															 
	return P -> PointInPolygon_wn(P, V; kwargs..., prepare_vertices=false)

end 








#===========================================================================#
#
# Generate vertices of a d-dimensional body based on d vectors
# (1D: segmemt, 3D: parallelipiped, 2D:parallelogram)
#
#---------------------------------------------------------------------------#


function rhomboidVertices_fromDVectors(v::AbstractMatrix; dim)::Matrix

	Utils.CombsOfVecs10(v, 1, 0; dim=dim)

end

# hypercube 

#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#



#function points_common_type(ps...)::DataType 
#
#	T = typeof(ps[1][1])  
#
#	for p in ps 
#
#		T = promote_type(T,typeof(p[1]),typeof(p[2]))
#		
#	end 
#
#	return GeometryBasics.Point2{T}
#
#end 
#
#function Points(V::AbstractMatrix; dim=2)::Vector{<:GeometryBasics.Point2}
#
#	Points(eachslice(V, dims=dim)...)
#
#end 
#
#function Points(ps::Vararg{<:AbstractVector}; kwargs...)::Vector{<:GeometryBasics.Point2} 
#
#	P = points_common_type(ps...) 
#
#	return P[p isa P ? p : P(p) for p in ps]
#
#end 
#
#
#
#function Points(ps::AbstractVector{<:GeometryBasics.Point2}; kwargs...)::Vector{<:GeometryBasics.Point2} 
#
#	P = points_common_type(ps) 
#
#	for p_ in ps 
#
#		p_ isa P && continue 
#
#		return Vector{P}[p isa P ? p : P(p) for p in ps]
#
#	end 
#
#	return ps
#
#end  
#









#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#



function concave_hull(points::AbstractVector{<:AbstractVector{T}}; outdim=2)::Vector{Vector{T}} where T<:Real
	
	ch = ConcaveHull.concave_hull(points)


	return ch.vertices 

#	s = ArrayOps.Array_from_ListIndsVals(
#							[outdim,[2,1][outdim]],
#							[length(ch.vertices), length(ch.vertices[1])]
#							)
#
#	@show ch 
#@show s  
#
#	out = similar(ch, Tuple(s))
#
#	
#	for (i,v) in enumerate(ch.vertices)
#
#		setindex!(selectdim(out, dim, i), v, :)
#
#	end 
#
#
#	return out

end 



function concave_hull(points::AbstractMatrix{T}; dim=2, kwargs...)::Matrix{T} where T<:Real

	ch = concave_hull(collect(eachslice(points, dims=dim)); kwargs...)

 	return cat(Utils.VecAsMat.(ch,dim)...,dims=dim)


end




#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function convex_hull_segm(vertices::AbstractVector{<:AbstractVector{<:Real}}
													)::Vector{LazySets.LineSegment}

	hull = LazySets.convex_hull(vertices) 

	return [LazySets.LineSegment(hull[i-1 + (i<2)*length(hull)],hull[i]) for i=axes(hull,1)]

end 




#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

#Fix2(f, x) behaves similarly to y->f(y, x).
#Fix1(f, x) behaves similarly to y->f(x, y).


function nearest_point(x::Union{<:AbstractVector{<:Real},
																<:LazySets.Line,
																},
											 points::AbstractVector{T}
											 )::T where T<:AbstractVector{<:Real}

# distance must use Line, not Line2D 
# project uses Line2D, not line 

	argmin(Base.Fix1(LazySets.distance, x), points)

end  

function project_point_to_line(p::AbstractVector{<:Real},
															 line::Union{LazySets.Line,LazySets.Line2D}, 
															 )::Vector{Float64} 
	
	LazySets.project(p, line_to_line2D(line))

end 

															 
function nearest_points(line::Union{LazySets.Line,LazySets.Line2D}, 
											 points::AbstractVector{T},
											 )::Tuple{AbstractVector{<:Real},T
																} where T<:AbstractVector{<:Real}

# distance must use Line
# project uses Line2D 

	p = nearest_point(line, points)

	return project_point_to_line(p, line), p

end 



function line_to_line2D(line::LazySets.Line)::LazySets.Line2D 

	LazySets.Line2D(line.p, line.p+line.d)

end    



function line_to_line2D(line::LazySets.Line2D)::LazySets.Line2D 

	line 

end  



function any_intersects(line::Union{LazySets.Line,
																		LazySets.Line2D,
																		LazySets.LineSegment},
												objects::AbstractVector)::Bool

#	!all(Base.Fix2(LazySets.is_intersection_empty,line),objects) 

	for obj in objects 

		LazySets.is_intersection_empty(line, obj) || return true 

	end 

	return false 
	
end 
#
#
#function any_intersects(objects::AbstractVector)::Function
#
#	Base.Fix2(any_intersects, objects)
#
#end 

function find_close_points(x::Union{<:AbstractVector{<:Real},
																<:LazySets.Line,
																},
											 points::AbstractVector{<:AbstractVector{<:Real}};
											 atol::Real=1e-8
											 )::Vector{Int}

	findall(<(atol) ∘ Base.Fix1(LazySets.distance, x), points)


end 





#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#



function Maximize_ContactSurface(starting_points::AbstractMatrix{<:Real}, 
																 direction::AbstractVector{<:Real}, 
																 vertices::AbstractMatrix{<:Real}; 
																 dim::Int,
																 kwargs...)::Vector{Float64}

	Maximize_ContactSurface(
										Vector{Float64}.(eachslice(starting_points,dims=dim)), 
										convert(AbstractVector{Float64},direction),
										collect(eachslice(vertices,dims=dim));
										kwargs...
										)

end 


function Maximize_ContactSurface(
								starting_points::AbstractVector{<:AbstractVector{<:Real}},
								direction::AbstractVector{Float64},
							 vertices::AbstractVector{<:AbstractVector{<:Real}};
																 kwargs...)::Vector{Float64}


	N = length(starting_points)

	lines = Vector{LazySets.Line}(undef,N)

	shifts = zeros(length(first(starting_points)),N)

	nrs = zeros(Int, N) 
	
	



	for i=1:N 


		lines[i] = LazySets.Line(starting_points[i], direction)


		
		shifts[:,i] .= nearest_point(lines[i], vertices)
		#closest vertex 


		shifts[:,i] -= project_point_to_line(selectdim(shifts,2,i), lines[i])
		# vector from the line to the closest vertex  

	end 

	
	hull_segm = convex_hull_segm(vertices) 


	for i=1:N,j=1:N 
		
		nrs[i] += any_intersects(LazySets.translate(lines[j], selectdim(shifts,2,i)),hull_segm)
			# count how many shifted lines intersect the hull 

	end 



	println.(eachcol(shifts))

	# choose the smallest shift that maximizes the number of rays inside 

	shift_out = copy(argmin(LA.norm, 
													(selectdim(shifts, 2, i) for i=findall(==(maximum(nrs)),nrs))))



	for i=1:N 

		# implement the shift so far 
		LazySets.translate!(lines[i], shift_out)
		
		shifts[:,i] .= 0 #-shift_out 
		shifts[:,i] .-= starting_points[i]


		# vertices extremely close to the shifted line 
		inds_ps = find_close_points(lines[i], vertices)

		if isempty(inds_ps)
			
			nrs[i] = 0  

		else 
		
			nrs[i] = 1 



			shifts[:,i] += nearest_point(-selectdim(shifts,2,i)+shift_out, 
																	 view(vertices, inds_ps)) 


				
		end 

	end   



	shift_out += argmin(LA.norm, 
											eachcol(selectdim(shifts, 2, findall(==(1),nrs)))
											)
	# smallest parallel shift 


#	shift_out -= direction # no contact 


	return shift_out 


end 


#function Maximize_ContactSurface2(starting_points::AbstractMatrix{<:Real}, 
#																 direction::AbstractVector{<:Real}, 
#																 vertices::AbstractMatrix{<:Real}; 
#																 dim::Int,
#																 kwargs...)::Vector{Float64}
#
#
#	dir_parallel = LinearAlgebra.normalize(direction) 
#
##	D = [dir_parallel;; dir_parallel[2]; -dir_parallel[1]] # on columns 
#
#	D = [dir_parallel[1]; -dir_parallel[2] ;;
#			 dir_parallel[2]; dir_parallel[1]
#			 ] # on rows 
#
#
#	m = Maximize_ContactSurface2(D*starting_points, D*vertices; kwargs...)
#
#		Maximize_ContactSurface2(D*transpose(starting_points),
#														 D*transpose(vertices);
#														 kwargs...)
#
#
#
#end 
































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































#############################################################################
end
