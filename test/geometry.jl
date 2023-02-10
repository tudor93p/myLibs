import myLibs: Geometry, Utils,Algebra,Lattices
import PyPlot; const plt=PyPlot
import LinearAlgebra; const LA=LinearAlgebra
import Random, LazySets


Random.seed!(8)

PyPlot.close() 

dim = 1 
dim2 = 2 

each(A) =eachslice(A,dims=dim)

scatter(A;kwargs...) = plt.scatter(eachslice(A,dims=dim2)...; kwargs...)

starting_points = Utils.PathConnect((a->vcat(a,a + [0 4]))(rand(1,2)), 5; dim=dim)[1]



direction = [1.0,0.0]#rand(2) #[1.0,1.0] 

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





	



scatter(starting_points,label="start")



scatter(V,label="vert")


	X = extrema(vcat(V,starting_points))
	plt.gca().autoscale(enable=false)


points = collect.(eachrow(V))

	rl = LazySets.Line(;from=rand(2).-0.5,to=rand(2).-0.5)

@testset "nearest point" begin 

	@test Geometry.nearest_point(points[5],points)≈points[5] 

	p1,p2 = p12 = rand(points,2)

#	scatter(hcat(p12...)',c="green")



	Y = [rl.p[2] + (x-rl.p[1])*rl.d[2]/rl.d[1] for x in X]


	p0 = Geometry.nearest_point(rl,points)


	@test 	p0≈points[	argmin(LA.norm.([p0] .- points))]


#	plt.plot(X,Y,c="green")
#	plt.scatter(p0...,c="k",marker="x")
#	scatter(hcat(Geometry.nearest_points(rl,p12)...)',c="green")

	for i in Geometry.find_close_points(rl,points,atol=2)

#		plt.scatter(points[i]...,c="r")

	end 



end 


function gety(line::LazySets.Line,x)

	line.p[2] + (x-line.p[1])*line.d[2]/line.d[1] 
	
end 

function pltline(line, x; kwargs...)
	
	plt.plot(x, [gety(line,x1) for x1 in x]; kwargs...)
	
end 


l1 = LazySets.Line(starting_points[1,:],direction)

pltline(l1, X,c="b")







m = Geometry.Maximize_ContactSurface(starting_points, direction, V; dim=dim)

																		 
#sp = vcat((hcat((p+m)...) for p in each(starting_points))...)
#sp = vcat((hcat((p+m+direction)...) for p in each(starting_points))...)

#while any(Geometry.PointInPolygon_wn(V; dim=dim), each(sp))
#
#	for i in axes(sp, dim)
#
#			selectdim(sp, dim, i) .-= direction*dmax 
#
#	end 
#
#end 
#

#scatter(sp;c="red",label="shifted")

@show m  



pltline(LazySets.Line(starting_points[1,:]+m,direction),X,c="red")
plt.gca().autoscale(enable=true)
	plt.scatter(eachcol(starting_points .+ transpose(m))...,c="red",marker="x")

plt.legend()
plt.show()


