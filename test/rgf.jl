import myLibs: GreensFunctions, Lattices, H_Superconductor,TBmodel, LayeredSystem,Algebra, ObservablesFromGF,Utils, BandStructure
import PyPlot, LinearAlgebra


# tests as in
# J Comput Electron (2013) 12:203–231 
# The recursive Green’s function method for graphene
# Caio H. Lewenkopf · Eduardo R. Mucciolo




function plot_atoms(n::Int, atoms::AbstractMatrix; kwargs...)

	PyPlot.figure(n)

	plot_atoms(PyPlot.gca(), atoms; kwargs...)

end  

function plot_atoms(ax, atoms::AbstractMatrix; 
										max_atoms::Int=9999,#40,
										kwargs...)
	
	nr_atoms =  Lattices.NrVecs(atoms)


	plot_at = Lattices.Vecs(atoms, 1:min(nr_atoms,max_atoms))

	ax.scatter(eachrow(plot_at)...; kwargs...) 

	ax.set_aspect(1)  


	if nr_atoms>max_atoms 

		ax.scatter(Lattices.Vecs(atoms, max_atoms+1)...;
							 Utils.dict_diff(kwargs,:c,:color)...,
							c="red")
	end 

end  

function plot_layers(ax_or_n; max_nr_layers::Int=999,#30, 
										 AtomsOfLayer::Function,
										 NrLayers::Int,
										 IndsAtomsOfLayer=nothing,
										 LayerOfAtom=nothing,
										 kwargs...)

	for layer in 1:NrLayers
	
		atoms = AtomsOfLayer(layer)
	
		per_layer = Lattices.NrVecs(atoms)
	
		kw = layer==1 ? Dict(:label=>"Device") : ()
	
		color = isodd(layer) ? "black" : "gray" 
	
		plot_atoms(ax_or_n, atoms;
							 color = layer>max_nr_layers ? "red" : color,
							 Utils.dict_diff(kwargs,:c,:color)...,
							 kw...)
	
		layer>max_nr_layers && break 

end 


end 





function honeycomb_x_armchair()

	latt = Lattices.HoneycombLattice() 
	
	Lattices.Superlattice!(latt, [[1 1]; [-1 1]])  

	return latt 

end 


function make_y_ribbon!(latt::Lattices.Lattice, NrAtomsPerLayer::Int)

	Lattices.Superlattice!(latt, [NrAtomsPerLayer,1])
	
	Lattices.KeepDim!(latt, 2, :absolute)

	return latt 
end 


function ribbon_x_armchair2(NrAtomsPerLayer::Int)

	L0 = honeycomb_x_armchair()
	
	#four_atoms[:,argmin(four_atoms[1,:])] .+= latt.LattVec[:,2]

	four_atoms = Lattices.sortvecs(Lattices.PosAtoms(L0);	
																 by=LinearAlgebra.norm
																 ) |> Lattices.eachvec


	@assert length(four_atoms)==4

	L = Lattices.Lattice(L0.LattVec,
											 [string(i)=>v for (i,v) in enumerate(four_atoms)],
											 nothing,
											 L0.LattDims
											 )

	return make_y_ribbon!(L, NrAtomsPerLayer)

end 


function slices_armchair(NrLayers::Int, NrAtomsPerLayer::Int) 

	L = ribbon_x_armchair2(NrAtomsPerLayer)


	function AtomsOfLayer(layer::Int)::Matrix{Float64}

		n = Int(floor((layer-1)/4))

		return Lattices.Atoms_ManyUCs(L; Ns=n, label=string(layer-4n))
	
	end 
	
	
	function IndsAtomsOfLayer(layer::Int)::Vector{Int}
	
		@assert 1<=layer<=NrLayers
	
		NrAtomsPerLayer*(layer-1) .+ (1:NrAtomsPerLayer)
	
	end 
	
	
	function LayerOfAtom(atom::Int)::Int 
	
		@assert  1<=atom<=NrLayers*NrAtomsPerLayer
	
		div(atom-1,NrAtomsPerLayer) + 1 
	
	end 
	
	
	for l in 1:NrLayers
		
		I = IndsAtomsOfLayer(l)
		
		@assert size(AtomsOfLayer(l),2)==length(I)
	
		for a in I 
	
		@assert LayerOfAtom(a)==l
	end 
	end  
	
	return Dict(

			:NrLayers=> NrLayers,

			:LayerOfAtom => LayerOfAtom,

			:IndsAtomsOfLayer => IndsAtomsOfLayer,

			:AtomsOfLayer => AtomsOfLayer

			)

end



function square_latt_lead(contact::Int, dir::Int, label::AbstractString;
													AtomsOfLayer::Function, 
													square_uc::Int=1,
													kwargs...)

	atoms = AtomsOfLayer(contact) 


	a = dir * (if Lattices.NrVecs(atoms)>1 
							 
							 only(Lattices.Distances(atoms,1))

							else 1

							end)

	atoms[1,:] .+= a

	latt = Lattices.SquareLattice(a, label, atoms, nothing, 1)

	square_uc==4 && Lattices.Superlattice!(latt, 4)
	
	return latt, abs(a)

end  


function honeycomb_latt_lead(contact::Int, dir::Int, label::AbstractString;
														 AtomsOfLayer::Function, 
														 kwargs...)

	atoms = mapreduce(AtomsOfLayer, Lattices.Cat,
										contact+dir :dir: contact+4dir)

	L = ribbon_x_armchair2(Lattices.NrVecs(AtomsOfLayer(1)))

	a = only(Lattices.Distances(atoms,1))

	return Lattices.Lattice(dir*L.LattVec,atoms,
													nothing,
													L.LattDims), a
end 

function lead_contacts(contact::Int,
											 lead::Lattices.Lattice, a::Real; 
											 AtomsOfLayer::Function, IndsAtomsOfLayer::Function, 
											 kwargs...)


	isBond = Algebra.EuclDistEquals(a; tol=1e-5, dim=2)

	lc = LayeredSystem.get_LeadContacts(AtomsOfLayer(contact), lead, isBond)

	return IndsAtomsOfLayer(contact)[lc]

end 

two_lead_args(;NrLayers::Int, kwargs...) =  [(1, -1, "Left"), (NrLayers, 1, "Right")] 



function device_atoms(;AtomsOfLayer::Function, NrLayers::Int, kwargs...)

	mapreduce(AtomsOfLayer, hcat, 1:NrLayers)

end 

dev_hopp = H_Superconductor.SC_Domain((ChemicalPotential=0,),	[1]) 

get_TBL(l::Lattices.Lattice) = Lattices.NearbyUCs(l, 1)


function prep_lead(label, latt, lead_hopp, del=delta)

	intra, inter = TBmodel.BlochHamilt_ParallPerp(latt, get_TBL, lead_hopp)
	
	function gf(E::Real; k=[])::Matrix{ComplexF64}
	
		GreensFunctions.GF_SanchoRubio(E + del, intra(k), only(values(inter)), "+")
	
	end
	
	
	lead_hopp_matr(args...) = TBmodel.HoppingMatrix(args...; lead_hopp...)

	return LayeredSystem.PrepareLead(label, latt, lead_hopp_matr, lead_hopp_matr, gf)
	
end 

function get_lead_hopp(latt) 
	
	H_Superconductor.SC_Domain((ChemicalPotential=0,),
														 Lattices.Distances(latt,1))

end


colors = [["orange","green"],["blue","gold"]]


nr_atoms_layer = 1+ 3*2

#nr_layers = Int(ceil(nr_atoms_layer*5.2))

#nr_layers = Int(ceil(nr_atoms_layer*2))
#

nr_layers = 120 


delta = 0.01im

@show nr_layers nr_atoms_layer

LayerAtomRels = slices_armchair(nr_layers, nr_atoms_layer) 






function set_title(fig,nr_layers=nr_layers,nr_atoms_layer=nr_atoms_layer)

	fig.suptitle(string("Figure ",
										 fig.number,
										 ": $nr_layers slices x $nr_atoms_layer atoms"))

end





#
#inds_bonds = NTuple{2,Int}[] 
#Rs_bonds = Vector{Vector{Float64}}[]
#
#for l1 in 1:NrLayers 
#
#	NrAtomsPerLayer > 50 && continue 
#
#	I1,A1 = IndsAtomsOfLayer(l1), AtomsOfLayer(l1)
#
#	for l2 in l1:min(l1+1,NrLayers)
#
#		I2,A2 = IndsAtomsOfLayer(l2), AtomsOfLayer(l2)
#
#		for (i1,i2) in Algebra.get_Bonds(A1,A2,1.0;dim=2,order_pairs=(l1==l2))
#
#			a1,a2 = A1[:,i1], A2[:,i2] 
#
#
#			PyPlot.plot([a1[1],a2[1]], [a1[2],a2[2]], c="k")
#
#			@assert isapprox(1, LinearAlgebra.norm(a1-a2))
#
#			push!(inds_bonds, (I1[i1],I2[i2]))
#			push!(Rs_bonds, [a1,a2])
#
#
#		end 
#
#	end 
#end 
#
#@show length(inds_bonds) 
#
#order = sortperm(inds_bonds)
#
#inds_bonds=inds_bonds[order]
#Rs_bonds=Rs_bonds[order]
#
#		
#
#
function NewGeometry(LayerAtomRels, LeadContacts, leads)

	NG = LayeredSystem.NewGeometry(rand(0,0), 
															 LayerAtomRels; 
															 LeadContacts=LeadContacts,
															 Leads=leads, dim=2, 
															 dev_hopp...)

	@assert length(NG)==4 
	
	return NG  

end


#LayerAtom_adv, Slicer_adv, LeadRels_adv, VirtLeads_adv = NG_adv

function GF((LayerAtom,Slicer,LeadRels,VirtLeads))

	GreensFunctions.GF_Decimation(dev_hopp, VirtLeads, Slicer;
																			 LayerAtom...)

end 

function get_SE(g, (LayerAtom,Slicer,LeadRels,VirtLeads))

	GreensFunctions.SelfEn_fromGDecim(g, VirtLeads, Slicer)

end 
function setIndRe!(X,x,inds...)

	@assert imag(x)<1e-10  x

	X[inds...] = real(x)

end 

#
#
##ga = Ga(Energy) 
#
#	SVTs = map([0,0.3]) do Energy
#	
#		gr1 = Gr(Energy)
#		
#		get_SE_ret1 = GreensFunctions.SelfEn_fromGDecim(gr1, VirtLeads_ret, Slicer_ret)
#		
#		
#		return ObservablesFromGF.SiteTransmission(gr1,
#																						 [ones(ComplexF64,1,1) for b in inds_bonds],
#																						 inds_bonds,
#																						 Rs_bonds,
#																						 get_SE_ret1,
#																						 "Left",2; dim=2)
#	end 
#	
#	@show size.(SVTs) 
#	
#	
#	
#	
#	XY = zero(SVTs[1]) 
#	
#	
#	for L in 1:NrLayers 
#	
#		for (i,r) in zip(IndsAtomsOfLayer(L), eachcol(AtomsOfLayer(L)))
#	
#			XY[:,i]=r
#	
#		end 
#	
#	end 
#	
#	
#	
#	fig,axes = PyPlot.subplots(1,2,num=8)
#	
#	for (svt,ax,col) in zip(SVTs,axes,["red","green"])
#	
#		ax.quiver(eachrow(XY)..., eachrow(svt)...,color=col,zorder=10)
#	
#		plot_atoms(ax,XY,color="k",zorder=0,s=10)
#	
#	end 
#	
#	
#
#
#end 
#
#
#
#
#
#Y = zeros(4, length(ENERGIES)) 
#
#labels = Vector{String}(undef,4)

################# 




#PyPlot.figure(6)
#
#y1 = Utils.Rescale(Y[1,:],[1,1.7],[0,maximum(Y[1,:])])
#
#PyPlot.plot(y1,ENERGIES,alpha=0.6,c="k",label="Caroli")
#PyPlot.plot(y1,-ENERGIES,alpha=0.6,c="k")
#PyPlot.legend()












