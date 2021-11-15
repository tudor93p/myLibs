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
										max_atoms::Int=100,
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

function plot_layers(ax_or_n; max_nr_layers::Int=30, 
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




function honeycomb_x(termination_x::Symbol)

	Lattices.HoneycombLattice_xy("x"; termination_x=termination_x)

end 



function make_y_ribbon!(latt::Lattices.Lattice, NrAtomsPerLayer::Int)

	Lattices.Superlattice!(latt, [1, NrAtomsPerLayer])
	
	Lattices.KeepDim!(latt, 1, :absolute)

	return latt  

end 


function layer_width(term_x::Symbol, N::Int)::Int 

	term_x==:zigzag && return Int(round(N/2))

	term_x==:armchair && return Int(3round((N-1)/3)+1)

end


function ribbon_x(term_x::Symbol, NrAtomsPerLayer::Int)

	make_y_ribbon!(honeycomb_x(term_x), 
								 layer_width(term_x, NrAtomsPerLayer))

end 




function slices_ribbon(term_x::Symbol, NrLayers::Int, NrAtomsPerLayer::Int)

	LayeredSystem.LayerAtomRels((ribbon_x(term_x, NrAtomsPerLayer),
																NrLayers), "sublatt";test=true)[1]

end











function square_latt_lead(contact::Int, dir::Int, label::AbstractString;
													AtomsOfLayer::Function, 
													four_uc::Bool=false,
													kwargs...)

	atoms = AtomsOfLayer(contact) 


	a = dir * (if Lattices.NrVecs(atoms)>1 
							 
							 only(Lattices.Distances(atoms,1))

							else 1

							end)

	atoms[1,:] .+= a

	latt = Lattices.SquareLattice(a, label, atoms, nothing, 1)

	four_uc && Lattices.Superlattice!(latt, 4)
	
	return latt, abs(a)

end  


function honeycomb_latt_lead(contact::Int, dir::Int, label::AbstractString;
														 term_x::Symbol,
														 AtomsOfLayer::Function, 
														 four_uc::Bool=false,
														 kwargs...)
	
	nr_uc = term_x==:armchair ? 4 : 2

	atoms = mapreduce(AtomsOfLayer, Lattices.Cat,
										contact+dir :dir: contact+nr_uc*dir)



	L = ribbon_x(term_x, Lattices.NrVecs(AtomsOfLayer(1)))

	a = only(Lattices.Distances(atoms,1))

	latt = Lattices.Lattice(dir*L.LattVec,atoms,
													nothing,
													L.LattDims)
	
	four_uc && term_x==:zigzag && Lattices.Superlattice!(latt, 2) 

	return latt,a


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

dev_hopp = H_Superconductor.SC_Domain((ChemicalPotential=0,
																			 Peierls=(0.0, 3sqrt(3)/2, :x),
																			 ), [1]) 

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


function get_Bonds(;NrLayers::Int, AtomsOfLayer::Function,
									 IndsAtomsOfLayer::Function,
										kwargs...)

	dev_at = device_atoms(;NrLayers=NrLayers, AtomsOfLayer=AtomsOfLayer,
												IndsAtomsOfLayer=IndsAtomsOfLayer,
												kwargs...)

#	N = min(NrLayers, Int(ceil(max_atoms/length(IndsAtomsOfLayer(1)))))


	inds_bonds = Algebra.get_Bonds(AtomsOfLayer, 1.0, 
																IndsAtomsOfLayer, 
																[(l1,l2) for l1 in 1:NrLayers for l2 in l1:min(l1+1,NrLayers)];
																dim=2)

	Rs_bonds = Algebra.bondRs_fromInds(inds_bonds, dev_at; dim=2)


	return inds_bonds, Rs_bonds

end 
function plot_bonds(ax, Rs_bonds; max_atoms::Int=500, kwargs...)

	ax.set_aspect(1)

	for (a1,a2) in Rs_bonds[1:min(max_atoms,end)]

		ax.plot([a1[1],a2[1]], [a1[2],a2[2]]; kwargs...)
	
	end 

end 



colors = [["orange","green","peru","dodgerblue"],["blue","gold","darkviolet","forestgreen"]]


nr_layers = 20

nr_atoms_layer = 3*4+1 

#@assert iseven(nr_atoms_layer)

#nr_layers = Int(ceil(nr_atoms_layer*5.2))

#nr_layers = Int(ceil(nr_atoms_layer*2))
#



delta = 0.005im

#@show nr_layers nr_atoms_layer

#LayerAtomRels_armchair = slices_ribbon(:armchair, nr_layers, nr_atoms_layer)  
#
#LayerAtomRels_zigzag = slices_ribbon(:zigzag, nr_layers, nr_atoms_layer)



function set_title(fig,nr_layers=nr_layers,nr_atoms_layer=nr_atoms_layer)

	fig.suptitle(string("Figure ",
										 fig.number,
										 ": $nr_layers slices x $nr_atoms_layer atoms"))

end



function xy_PHS( E::AbstractVector{<:Real},
									f::AbstractVector{<:Real},
									E_on_axis::Int=1)


	E2 = vcat(-reverse(E),E)
	f2 = vcat(reverse(f),f)

	E_on_axis==1 && return E2,f2  

	E_on_axis==2 && return f2,E2

error()





end 






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

	@assert imag(x)<1e-9  x

	X[inds...] = real(x)

end 

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












