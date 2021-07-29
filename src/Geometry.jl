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



function Maximize_ContactSurface(starting_points::AbstractMatrix{Float64}, 
																 direction::AbstractVector{Float64}, 
																 vertices::AbstractMatrix{Float64}; 
																 dim::Int,
																 prepare_vertices::Bool=true, 
																 kwargs...)::Vector{Float64}

	prepare_vertices && return Maximize_ContactSurface(starting_points, direction, prepare_polygon_vertices(vertices; dim=dim, kwargs...); dim=dim, prepare_vertices=false, kwargs...)


	poly_verts = Points(vertices; dim=dim)

	poly_sides = Segments(CGAL.convex_hull_2(poly_verts)) 



	shifts = map(eachslice(starting_points, dims=dim)) do start 

		line = Line(Points(start, start + direction)...)

#		any(side->do_intersect(line, side), poly_sides) && return zeros(2)

		p = partialsort(poly_verts, 1, by = p->squared_distance(line, p))

		return xyz(p - projection(line, p)) # from the line to the vertex  

	end 

	shift_perp = partialsort!(shifts, 1, by=LA.norm, rev=true) 
					# largest shift => all starting points happy  


	shifts = map(eachslice(starting_points, dims=dim)) do start_

		start = start_ + shift_perp 

		pstart,pstop = Points(start, start+direction)

		line = Line(pstart, pstop)

		ps = filter(p->squared_distance(line, p)<1e-12, poly_verts)

		isempty(ps) && return nothing 

		p = partialsort(ps, 1, by = p->squared_distance(pstart, p))

		return xyz(p-pstart) - direction

	end 

	shift_parallel = partialsort!(filter!(!isnothing,shifts), 1, by=LA.norm)


	return shift_perp + shift_parallel

end 




#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function ClosestContact_LinesPolygon(starting_points::AbstractMatrix{Float64},
																		 direction::AbstractVector{Float64},
																		 vertices::AbstractMatrix{Float64}; 
																		 dim::Int,
																		 prepare_vertices::Bool=true,
																		 kwargs...)::Vector{Float64}



	prepare_vertices && return ClosestContact_LinesPolygon(starting_points, direction, prepare_polygon_vertices(vertices; dim=dim, kwargs...); dim=dim, prepare_vertices=false, kwargs...)



	poly_verts = Points(vertices; dim=dim)

	
	lines = [Line(Points(start, start + direction + 1e-12rand(2))...) for start in eachslice(starting_points, dims=dim)]


	D = Algebra.OuterBinary(lines, poly_verts, squared_distance)  


	@show minimum(D)



#
#	argmin(eachrow(D))


#	x = map(eachslice(starting_points, dims=dim)) do start 
#
#		line = Line(Points(start, start + direction + 1e-12rand(2))...)
#
#		any(side->do_intersect(line, side), poly_sides) && return zeros(2)
#
#		p = partialsort(poly_verts, 1, by = p->squared_distance(line, p))
#
#		return xyz(p - projection(line, p)) # from the line to the vertex 
#
#	end 
#
#	println.(x)

	return [1.0]




end 
































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































#############################################################################
end
