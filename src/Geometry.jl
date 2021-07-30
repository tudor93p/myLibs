module Geometry
#############################################################################

import ConcaveHull, CGAL
#import GeometryBasics, 

import LinearAlgebra; const LA=LinearAlgebra

#using GeometryBasics: Point, MultiPoint, LineString, Line 
import GeometryBasics 

using CGAL: Point2 as Point, Line2 as Line, Segment2 as Segment, do_intersect, squared_distance, projection

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



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function Points(V::AbstractMatrix; dim=2)::Vector{Point}

	Points(eachslice(V, dims=dim)...)

end 


function Points(ps::Vararg{<:AbstractVector}; kwargs...)::Vector{Point}

	[Point(p[1],p[2]) for p in ps]

end 

function Points(ps::AbstractVector{<:Point}; kwargs...)::Vector{Point} 

	ps 

end 

function Points(ps::Vararg{<:Point})::Vector{Point}

	collect(ps)

end 


function Segments(ps::Vararg{<:Point,N})::Vector{Segment} where N

	[Segment(ps[i],ps[i+1]) for i in 1:N-1]

end  

function Segments(ps...; dim=2)::Vector{Segment}

	Segments(Points(ps...; dim=dim)...)
	

end 


function xyz(P)::Vector{Float64}

	[getproperty(CGAL, f)(P) for f in [:x, :y, :z][1:CGAL.dimension(P)]]

end

		

#function SquaredDistance(p1::Point, p2::Point)::Float64

#CGAL.Point2







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

function nearest_point(x::Union{<:Line, <:Point}, 
											 points::AbstractVector{Point})::Point 

	partialsort(points, 1, by = p->squared_distance(x, p)) 

end 

function nearest_points(line::Line, points::AbstractVector{Point})::Tuple{Point,Point}

	p = nearest_point(line, points)

	return projection(line, p), p

end 



function any_intersects(line::Line, objects::AbstractVector)::Bool

	any(x -> do_intersect(line, x), objects)

end 


function any_intersects(objects::AbstractVector)::Function

	anyint(line::Line)::Bool = any_intersects(line, objects)

end 




#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#



function Maximize_ContactSurface(starting_points::AbstractMatrix{Float64}, 
																 direction::AbstractVector{Float64}, 
																 vertices::AbstractMatrix{Float64}; 
																 dim::Int,
																 prepare_vertices::Bool=true, 
																 kwargs...)::Vector{Float64}

	prepare_vertices && return Maximize_ContactSurface(starting_points, direction, prepare_polygon_vertices(vertices; dim=dim, kwargs...); dim=dim, prepare_vertices=false, kwargs...)


	poly_verts = Points(vertices; dim=dim)



	line(start::AbstractVector) = Line(Points(start, start + direction)...)
	line(start::AbstractVector, shift::AbstractVector) = line(start+shift)

	
	iter_start = eachslice(starting_points, dims=dim) 


	ray_cuts = any_intersects(Segments(CGAL.convex_hull_2(poly_verts))) 

	nr_inside(args...)::Int = count((ray_cuts(line(s, args...)) for s in iter_start))


	shifts,nrs = Utils.zipmap(iter_start) do start

		p_line, p = nearest_points(line(start), poly_verts) 
								# the closest p and its projection on the line

		shift = xyz(p - p_line) 	# vector from the line to the vertex  

		return shift, nr_inside(shift)

	end 


	# choose the smallest shift that maximizes the number of rays inside 

	shift_perp = argmin(LA.norm, shifts[collect(nrs).==maximum(nrs)])



	shifts_parallel = Utils.mapif(!isnothing, iter_start .+ [shift_perp]) do start 

		L = line(start) 

		ps = filter(p->squared_distance(L, p)<1e-12, poly_verts)

		return isempty(ps) ? nothing : xyz(nearest_point(Point(start...), ps)) - start


	end 

									# smallest parallel shift 
	
	return shift_perp + argmin(LA.norm, shifts_parallel) - direction

end 


































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































#############################################################################
end
