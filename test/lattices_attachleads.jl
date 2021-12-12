import myLibs: Lattices, Geometry 

import PyCall 
using OrderedCollections:OrderedDict


println() 
	
path = "/media/tudor/Tudor/Work/scripts/python_libraries/tests/"

pushfirst!(PyCall.PyVector(PyCall.pyimport("sys")."path"), path)

pytest = PyCall.pyimport("Lattices_AlignLead_test")


println() 



function jl2(N::AbstractVector{<:Int}, 
						 polygon::AbstractMatrix{Float64}, 
						 desired_angle::Float64)

    println("jl2")

		SqLatt = Lattices.SquareLattice()
		
		Latt0 = Lattices.ShiftAtoms(Lattices.SquareLattice("Device"), n = (1 .-N)/2)

		Lattices.Superlattice!(Latt0, N) 
		
		Lattices.ReduceDim!(Latt0) 

		LattAtoms0 = Lattices.PosAtoms(Latt0)
		


		rmv_inds = map(!Geometry.PointInPolygon_wn(polygon, dim=1; order_vertices=true), eachcol(LattAtoms0))


		if all(rmv_inds) 

			rmv_inds[1:5] .= false

		end 

		d_nn = Lattices.Distances(SqLatt)[1]

		Latt = Lattices.RemoveAtoms(Latt0, rmv_inds, "A") 

		Lattices.RemoveSingleBonds!(Latt, d_nn)

		LattAtoms = Lattices.PosAtoms(Latt)

#		@show size(Lattices.SurfaceAtoms(LattAtoms, d_nn))

#		@show size(Geometry.concave_hull(LattAtoms; dim=2))

#		BA = Geometry.concave_hull(LattAtoms; dim=2)
    
		BA = Lattices.SurfaceAtoms(LattAtoms, d_nn)
    
    
	circle_outside, lead_pos = pytest.py1(desired_angle, transpose(LattAtoms0))
    
	Lead = Lattices.ShiftAtoms(Lattices.SquareLattice("Lead-initial"), 
														 lead_pos)
	


		LeadAtoms = Lattices.PosAtoms(Lead)


		return transpose(LattAtoms0), transpose(LeadAtoms), transpose(LattAtoms), Lead, SqLatt, transpose(BA)


end 

function jl3(Lead::Lattices.Lattice, 
						 sc::AbstractVecOrMat, 
						 LattAtoms::AbstractMatrix{Float64},
						 SqLatt::Lattices.Lattice)

	println("jl3")


	N = vcat(sc...)

	if length(N)==2

		Lattices.Superlattice!(Lead, [Int(i) for i in N])
		
	elseif length(N)==4 

		Lattices.Superlattice!(Lead, reshape(N, (2,2)))

	end 
			
	Lattices.ReduceDim!(Lead, 2) # python 1 = julia 2

	LeadAtoms2 = transpose(Lattices.PosAtoms(Lead))

	Lattices.Align_toAtoms!(Lead, transpose(LattAtoms))

	LeadAtoms3 = transpose(Lattices.PosAtoms(Lead))

	return LeadAtoms2, LeadAtoms3, zeros(0,2)
end  



Lead = Lattices.Lattice([-1.0 0.0; 0.0 2.0],
												OrderedDict{Any, Matrix{Float64}}("Lead-1" => [-4.5 -4.5; 1.5 2.5]),
												nothing,
												
												[1])
															 
Atoms = [-4.5 -4.5 -4.5 -4.5 -4.5 -3.5 -3.5 -2.5 -2.5 -1.5 -1.5 -0.5 -0.5 0.5 0.5 1.5 1.5 2.5 2.5 3.5 3.5 4.5 4.5 4.5 4.5 4.5; -2.0 -1.0 0.0 1.0 2.0 -2.0 2.0 -2.0 2.0 -2.0 2.0 -2.0 2.0 -2.0 2.0 -2.0 2.0 -2.0 2.0 -2.0 2.0 -2.0 -1.0 0.0 1.0 2.0];
															 
@show Lattices.PosAtoms(Lead)
@show Lattices.LattVec(Lead) 

@show size(Atoms) 


Lattices.Align_toAtoms!(Lead, Atoms)


@show Lattices.PosAtoms(Lead)

jl2([5,6],
		mapreduce(a->15*[cos(a) sin(a)], vcat,rand(5)*2pi),
		rand())


jl3(Lattices.ShiftAtoms(Lattices.SquareLattice(), -rand(2)),
		[[1 2];[1 0]],
		hcat(vcat.(Base.product(1:5,1:5)...)...) .- 1.0,
		Lattices.SquareLattice()
		)

#import PyPlot; fig,ax = PyPlot.subplots()
#
#pytest.plot0(14, 0* pi/180, 5, ax, jl2, jl3)
#

#pytest.plot(jl2, jl3) 







