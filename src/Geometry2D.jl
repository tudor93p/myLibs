module Geometry2D
#############################################################################

import ConcaveHull


import ..Utils





#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#



function concave_hull(points::AbstractVector{AbstractVector{<:Real}}; outdim=2)
	
	ch = ConcaveHull.concave_hull(points)


	s = Utils.Array_from_ListIndsVals(
							[outdim,[2,1][outdim]],
							[length(ch.vertices), length(ch.vertices[1])]
							)

	out = similar(ch, Tuple(s))

	
	for (i,v) in enumerate(ch.vertices)

		setindex!(selectdim(out, dim, i), v, :)

	end 


	return out

end 



function concave_hull(points::AbstractMatrix{<:Real}; dim=2, kwargs...)

	concave_hull(collect(eachslice(points, dims=dim)); kwargs...)

end

























































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































#############################################################################
end
