import myLibs: GreensFunctions, Lattices, H_Superconductor,TBmodel, LayeredSystem,Algebra, ObservablesFromGF,Utils, BandStructure
import PyPlot, LinearAlgebra

PyPlot.close.(1:10)

# tests as in
# J Comput Electron (2013) 12:203–231 
# The recursive Green’s function method for graphene
# Caio H. Lewenkopf · Eduardo R. Mucciolo


# Armchair # edges, M = 360 and N = 70 (aspect ratio W/L = 5.2) 
 

function plot_atoms(n::Int, atoms::AbstractMatrix; kwargs...)

	PyPlot.figure(n)

	plot_atoms(PyPlot.gca(), atoms; kwargs...)

end  

function plot_atoms(ax, atoms::AbstractMatrix; 
										max_atoms::Int = 30,
										kwargs...)
	
	nr_atoms =  Lattices.NrVecs(atoms)


	plot_at = Lattices.Vecs(atoms, 1:min(nr_atoms,max_atoms))

	ax.scatter(eachrow(plot_at)...; kwargs...) 

	ax.set_aspect(1)  


	if nr_atoms>max_atoms 

		ax.scatter(Lattices.Vecs(atoms, max_atoms+1)...;
							 kwargs..., color="red")

	end 

end  

function plot_layers(ax_or_n; max_nr_layers::Int=20, 
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
							 kw..., kwargs...)
	
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





#
#NrLayers =  360#N
#
#NrAtomsPerLayer = 1 + 3*23 # M 
#
#delta = 0.005im 
#
#
#default_LayerAtomRels = Dict(
#
#			:NrLayers=> NrLayers,
#
#			:LayerOfAtom => LayerOfAtom,
#
#			:IndsAtomsOfLayer => IndsAtomsOfLayer,
#
#			:AtomsOfLayer => AtomsOfLayer
#
#			)
#
#
#
#
#
#
##PyPlot.figure(7)
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
#leadlabels = ["Left","Right"]
#
#
#
#
#LeadLatts,LeadContacts = honeycomb_LLC 
#
#@show length.(LeadContacts)
#
#
#lead_hopp = H_Superconductor.SC_Domain((ChemicalPotential=0,),
#																			 Lattices.Distances(rand(LeadLatts),1))
#
#@show lead_hopp[:Hopping]([1,0],[0,0])
#@show lead_hopp[:Hopping]([0,0],[0,1]) 
#
#
#
#
#println()
#
#@show lead_hopp[:Hopping]([sqrt(3),0],[0,0])
#@show lead_hopp[:Hopping]([0,0],[0,sqrt(3)])
#
#println() 
#
#
#
#
#
#
##===========================================================================#
##
##
##
##---------------------------------------------------------------------------#
#
#
#
#
#
#
#
#
#
#dev_hopp = H_Superconductor.SC_Domain((ChemicalPotential=0,),	[1]) 
#		
#
#isBond = Algebra.EuclDistEquals(1; tol=1e-5,dim=2)
#
#
#
#
#
#
#lead_hopp_matr(args...) = TBmodel.HoppingMatrix(args...; lead_hopp...)
#
##@show lead_hopp_matr(AtomsOfLayer(1),Lattices.PosAtoms(LeadLatts[1]))
##@show lead_hopp_matr(AtomsOfLayer(NrLayers),Lattices.PosAtoms(LeadLatts[2]))
#
#get_TBL(l::Lattices.Lattice) = Lattices.NearbyUCs(l, 1)
#
#leads_ret,leads_adv = map([delta,-delta]) do del
#
#	map(zip(leadlabels, LeadLatts)) do (label, latt)
#
#		intra, inter = TBmodel.BlochHamilt_ParallPerp(latt, get_TBL, lead_hopp)
#	
#		function gf(E::Real; k=[])::Matrix{ComplexF64}
#	
#			GreensFunctions.GF_SanchoRubio(E + del, intra(k), only(values(inter)), "+")
#		
#		end
#	
#	
#		return LayeredSystem.PrepareLead(label, latt, lead_hopp_matr, lead_hopp_matr, gf)
#	
#	end 
#
#end
#
#
#
#println()
#
#
##===========================================================================#
##
## Lead lattice, spectrum and DOS 
##
##---------------------------------------------------------------------------#
#
#for l in LeadLatts 
#
#	k  = only(keys(l.Atoms))
#
#	v = Lattices.Atoms_ManyUCs(l,StopStart=(0,0))
#
#	plot_atoms(3, v; label=k, alpha=0.3, s=70)
#
#end 
#
#if false 
#
#	PyPlot.figure(6) 
#	
#	lead_Energies = LinRange(-4.2,4.2, 100) 
#	
#	lead_gf_labels = ["+","bulk"] 
#	
#	DOS = zeros(2, 2,length(lead_Energies))
#
#	dev_TBi = ([0],[[0,0]],mapreduce(AtomsOfLayer, hcat, 1:NrLayers))
#
#	spectrum_device = BandStructure.Diagonalize(TBmodel.Bloch_Hamilt(dev_TBi; dev_hopp...),zeros(0,1);dim=2)
#
#	PyPlot.scatter(spectrum_device["kLabels"], spectrum_device["Energy"], alpha=0.5, label= "Device",s=10,zorder=4)
#
#
#
#	for (il,(l1,L1)) in enumerate(zip(leads_ret,LeadLatts))
#	
#		H = TBmodel.Bloch_Hamilt(Lattices.NearbyUCs(L1,1);
#														 argH="k",
#														 lead_hopp...)
#	
#		kPoints, kTicks = Utils.PathConnect(Lattices.BrillouinZone(L1), 150; dim=2)
#	
#		spectrum = BandStructure.Diagonalize(H, kPoints; dim=2)
#	
#		PyPlot.scatter(spectrum["kLabels"], spectrum["Energy"], alpha=0.5,
#									 label= only(Lattices.sublatt_labels(L1)),s=10)
#	
#	
#		PyPlot.xlabel(spectrum["kLabel"])
#		PyPlot.ylabel("Energy")
#	
#		intra, inter = TBmodel.BlochHamilt_ParallPerp(L1, get_TBL, lead_hopp)
#	
#		intra_,inter_ = intra(), only(values(inter))
#	
#	
#	#	il==1 && @show intra_ inter_ 
#	
##		for e0 in LinearAlgebra.eigvals(Matrix(intra_))
#	
#	#		PyPlot.plot([0,1.7],[e0,e0],linestyle=":",c="gray")
#	
##		end 
#
#
#
##		@time begin 
#		for (ie,E) in enumerate(lead_Energies)
#		
#			gfs = GreensFunctions.GF_SanchoRubio(E + delta, intra_, inter_; 
#																					 target=join(lead_gf_labels))
#
##			@assert isapprox(GreensFunctions.GF_SanchoRubio(E + delta, intra_, inter_, "-"),
##											 GreensFunctions.GF_SanchoRubio_rec(E + delta, intra_, inter_;																			 target="-")["-"],
##											 rtol=1e-8)
##
##@assert isapprox(gfs["+"], l1[:GF](E)[1], rtol=1e-8)
#
#
#			for (igf,gflab) in enumerate(lead_gf_labels)
#
#				DOS[il,igf,ie] = ObservablesFromGF.DOS(gfs[gflab])
#		
#			end 
#		
#	
#			print("\r",round(100ie/length(lead_Energies),digits=1),"%     ") 
#
##		end 
#		end 
#
#		println()#"\r                     ")
#
#	
#	end 
#
#
#	DOS = Utils.Rescale(DOS, [1,1.7], [0, maximum(DOS)])  
#	
#	
#		
#	for (igf,gflab) in enumerate(lead_gf_labels)
#	
#		dos0 = DOS[1, igf, :]
#	
#		if all((isapprox(DOS[il, igf, :],dos0) for il in 2:length(leads_ret)))
#	
#			PyPlot.plot(dos0, lead_Energies; label=gflab, alpha=0.8)
#	
#		else 
#	
#			for (il,l1) in enumerate(leads_ret)
#			
#				PyPlot.plot(DOS[il, igf, :], lead_Energies;
#										label=string(l1[:label]," ",gflab), alpha=0.5)
#	
#			end 
#	
#		end 
#	
#	end 
#	
#	PyPlot.xlim(0,1.7)
#	PyPlot.legend() 
#	
#end 
##===========================================================================#
##
##
##
##---------------------------------------------------------------------------#
#
#														
#
#
#NG_ret = LayeredSystem.NewGeometry(rand(0,0), 
#															 default_LayerAtomRels; 
#															 LeadContacts=LeadContacts,
#															 Leads=leads_ret, dim=2, 
#															 dev_hopp...)
#
#NG_adv = LayeredSystem.NewGeometry(rand(0,0), 
#															 default_LayerAtomRels; 
#															 LeadContacts=LeadContacts,
#															 Leads=leads_adv, dim=2, 
#															 dev_hopp...)
#
#@assert length(NG_ret)==length(NG_adv)==4 
#
#
#LayerAtom_ret, Slicer_ret, LeadRels_ret, VirtLeads_ret = NG_ret
#LayerAtom_adv, Slicer_adv, LeadRels_adv, VirtLeads_adv = NG_adv
#
#
#Gr = GreensFunctions.GF_Decimation(dev_hopp, VirtLeads_ret, Slicer_ret;
#																			 LayerAtom_ret...)
#
#Ga = GreensFunctions.GF_Decimation(dev_hopp, VirtLeads_adv, Slicer_adv;
#																			 LayerAtom_adv...)
#
#
#if false 
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












