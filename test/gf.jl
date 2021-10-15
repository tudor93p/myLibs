import myLibs: GreensFunctions, LayeredSystem, Lattices, Algebra, Utils

P1 = (length = 60, Barrier_height = 1.75, Barrier_width = 0.03, SCDW_phasediff = 0.78, SCDW_p = 2, SCDW_width = 0.005, SCDW_position = 0, SCpx_magnitude = 0.4, delta = 0.002, width = 30, SCpy_magnitude = 0.4)

P2 = (Attached_Leads = "AB", Lead_Coupling = 1.0, Lead_Width = 0.5)

P3 = [(ChemPot = 0.0, Hopping = 1.0, Label = "A", Direction = -1, Contact = (-1, 1), Lead_Width = 15, Lead_Coupling = 1.0, SCbasis = true), (ChemPot = 0.0, Hopping = 1.0, Label = "B", Direction = 1, Contact = (1, -1), Lead_Width = 15, Lead_Coupling = 1.0, SCbasis = true)]




Dev = Lattices.SquareLattice()

NxNy = [P1[:length], P1[:width]]

Lattices.ShiftAtoms!(Dev, n=(1 .- NxNy)/2)

Lattices.Superlattice!(Dev, NxNy) 

Lattices.ReduceDim!(Dev)


DevAtoms = Lattices.PosAtoms(Dev) 

@show Lattices.NrVecs(DevAtoms)

SurfaceAtoms = Lattices.SurfaceAtoms(DevAtoms,1)

SurfaceAtoms2 = Lattices.filterAtoms_byNrBonds(-5:-1, DevAtoms, 1)


@assert isapprox(SurfaceAtoms,SurfaceAtoms2) 


@show Lattices.NrVecs(SurfaceAtoms)  

CornerAtoms = Lattices.filterAtoms_byNrBonds(-2, DevAtoms, 1)


@show Lattices.NrVecs(CornerAtoms)  

qs = Utils.Quadrant(CornerAtoms, Algebra.Mean(CornerAtoms, 2), dim=2) 

#@show qs 

CornerAtoms = CornerAtoms[:,sortperm(qs)]

#println.(eachcol(CornerAtoms))

LeadLatts = map(P3) do lead_params 

	d = lead_params[:Direction]::Int 

	sDir,Dir = sign(d),abs(d)

	L = Lattices.SquareLattice(string("Lead", sDir>0 ? "+" : "-", Dir))

	Lattices.KeepDim!(L, Dir) 


	Lattices.Superlattice!(L, lead_params[:Lead_Width])

	@show Lattices.NrVecs(Lattices.PosAtoms(L))

	return L

end 


	

#	!get_full && return L.Reduce_Dim(0,complement=true)
#	return L.Reduce_Dim(0,complement=true), (L,[1])




















#h = LayeredSystem.get_hoppings(Hopping[:Hopping], Slicer, VirtLeads)
#
#BondHoppings = [h(iR...)' for iR in zip(Bonds...)]
#


#G = GreensFunctions.GF_Decimation(
#				Hopping, VirtLeads, Slicer;
#				LayerAtom...)

