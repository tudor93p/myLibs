import myLibs: Geometry, Utils,Algebra,Lattices
import PyPlot; const plt=PyPlot
import LinearAlgebra; const LA=LinearAlgebra

dim = 1 
dim2 = 2 

each(A) =eachslice(A,dims=dim)

scatter(A;kwargs...) = plt.scatter(eachslice(A,dims=dim2)...; kwargs...)

starting_points = Utils.PathConnect((a->vcat(a,a + [0 4]))(rand(1,2)), 5; dim=dim)[1]



direction = [1.0,0.0] 

center = Algebra.Mean(starting_points, dim) + [0,5]

for i in axes(starting_points,dim)

	selectdim(starting_points,dim,i) .-= center 

end 








angles = sort(rand(8))*2pi 

#angles = range(0,2pi,length=21)
#
polygon = hcat(cos.(angles), sin.(angles))*9



V = vcat((hcat(a...) for a in Base.product(-5:5,-5:5))...) .+0.0

@show size(V) size(polygon)

@show size polygon


#plt.plot(eachslice(polygon, dims=dim2)...)

V = selectdim(V, dim, map(Geometry.PointInPolygon_wn(polygon, dim=dim),
												 each(V)))

scatter(V; c="lightgray")

V = transpose(Lattices.SurfaceAtoms(transpose(V),1))

#v,s = Algebra.Mean.([V, starting_points], dim) 



dmax = maximum(Algebra.OuterDist(V, V,dim=dim))



while any(Geometry.PointInPolygon_wn(V; dim=dim), each(starting_points))

	for i in axes(starting_points, dim)

			selectdim(starting_points, dim, i) .-= direction*dmax 

	end 

end 





	



scatter(starting_points)







scatter(V)











#plt.show()



m = Geometry.Maximize_ContactSurface(starting_points, direction*3dmax, V; dim=dim)
																		 
sp = vcat((hcat((p+m+direction)...) for p in each(starting_points))...)

while any(Geometry.PointInPolygon_wn(V; dim=dim), each(sp))

	for i in axes(sp, dim)

			selectdim(sp, dim, i) .-= direction*dmax 

	end 

end 


scatter(sp;c="red")



plt.show()


