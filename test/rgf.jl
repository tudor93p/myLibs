import myLibs: GreensFunctions, Lattices, H_Superconductor,TBmodel, LayeredSystem,Algebra, ObservablesFromGF,Utils, BandStructure
import PyPlot, LinearAlgebra

PyPlot.close.(1:10)

# tests as in
# J Comput Electron (2013) 12:203–231 
# The recursive Green’s function method for graphene
# Caio H. Lewenkopf · Eduardo R. Mucciolo

function plot_atoms(n::Int, atoms::AbstractMatrix; kwargs...)

	PyPlot.figure(n)

	PyPlot.gca().set_aspect(1) 

	PyPlot.scatter(eachrow(atoms)...; kwargs...)

end 

# Armchair # edges, M = 360 and N = 70 (aspect ratio W/L = 5.2) 
 
NrLayers = 50

NrAtomsPerLayer = 1



latt = Lattices.HoneycombLattice() 



for n in 0:3, m in 0:3

	#plot_atoms(2, Lattices.Atoms_ManyUCs(latt; Ns=[n,m]))

end 

for arrow in Lattices.eachvec(latt.LattVec) 

	#PyPlot.arrow(0,0,arrow...)
	
end   







Lattices.Superlattice!(latt, [[1 1]; [-1 1]]) 

# latt.Atoms["A"][:,1] .+= latt.LattVec[:,2]








#layer = unit_cell*nr_subs + sublattice 





atoms(n,m) = Lattices.Atoms_ManyUCs(latt; Ns=[n,m])


for n in 1:4

#	continue 

	#plot_atoms(1,hcat(atoms(n,0),atoms(0,n)))

end 

for arrow in Lattices.eachvec(latt.LattVec) 

#	PyPlot.arrow(0,0,arrow...)
	
end 


#PyPlot.scatter(eachrow(atoms(0,0))...,c="k")

for n in 0:2,m in 0:2
	
#	PyPlot.scatter(eachrow(atoms(m,n))...)
	
end 







four_atoms = Lattices.PosAtoms(latt) 


four_atoms[:,argmin(four_atoms[1,:])] .+= latt.LattVec[:,2]







L = Lattices.Lattice(latt.LattVec,
										 [string(i)=>v for (i,v) in enumerate(Lattices.eachvec(Lattices.sortvecs(four_atoms; by=first)))]
											
										 )

for (k,at) in L.Atoms 

#	PyPlot.scatter(eachrow(at)...,label=k)

end 
#
#mklab(n::AbstractVector) = string(n[2])



Lattices.Superlattice!(L, [NrAtomsPerLayer,1])#; Labels=mklab)#identity) 

Lattices.KeepDim!(L, 2)  



for n in 0:2

#	PyPlot.scatter(eachrow(Lattices.Atoms_ManyUCs(L; Ns=n))..., label=n)

end 




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

default_LayerAtomRels = Dict(

			:NrLayers=> NrLayers,

			:LayerOfAtom => LayerOfAtom,

			:IndsAtomsOfLayer => IndsAtomsOfLayer,

			:AtomsOfLayer => AtomsOfLayer

			)




for layer in 0:NrLayers+1 

	plot_atoms(3, AtomsOfLayer(layer); c=layer in 1:NrLayers ? "k" : "red")

end 








leadlabels = ["Left","Right"]



LeadLatts,LeadContacts = Utils.zipmap(zip([-1,1],leadlabels,[1,NrLayers])
																		 ) do (dir,label,contact)

	layer = AtomsOfLayer(contact) 



	ll = Lattices.Lattice(dir*L.LattVec,
												label=>mapreduce(AtomsOfLayer, hcat,
																				 contact+dir :dir: contact+4dir),
												nothing,
												L.LattDims)


	a = only(Lattices.Distances(ll,1))


#
#	atoms = AtomsOfLayer(contact + dir)
#
#	a = dir* (if size(atoms,2)>1 
#
#							only(Lattices.Distances(atoms,1))
#
#						else 
#
#							abs(layer[1,1]-atoms[1,1])
#
#						end)
#
#	atoms[1,:] .+= -atoms[1,1] + layer[1,1] + a 
#
#	ll = Lattices.SquareLattice(a, label, atoms, nothing, 1)
#
#	Lattices.Superlattice!(ll,4)





	isBond_lead = Algebra.EuclDistEquals(abs(a); tol=1e-5, dim=2)



	lc = LayeredSystem.get_LeadContacts(layer, isBond=isBond_lead,
														Leads=[Dict(:head=>[Lattices.PosAtoms(ll)])])

	return ll, IndsAtomsOfLayer(contact)[lc[1]]

end  .|> collect

@show LeadContacts 


lead_hopp = H_Superconductor.SC_Domain((ChemicalPotential=0,),
																			 Lattices.Distances(rand(LeadLatts),1))

@show lead_hopp[:Hopping]([1,0],[0,0])
@show lead_hopp[:Hopping]([0,0],[0,1]) 




println()

@show lead_hopp[:Hopping]([sqrt(3),0],[0,0])
@show lead_hopp[:Hopping]([0,0],[0,sqrt(3)])

println() 





error() 

#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#









dev_hopp = H_Superconductor.SC_Domain((ChemicalPotential=0,),	[1]) 
		

isBond = Algebra.EuclDistEquals(1; tol=1e-5,dim=2)


delta = 0.002im 




lead_hopp_matr(args...) = TBmodel.HoppingMatrix(args...; lead_hopp...)

@show lead_hopp_matr(AtomsOfLayer(1),Lattices.PosAtoms(LeadLatts[1]))
@show lead_hopp_matr(AtomsOfLayer(NrLayers),Lattices.PosAtoms(LeadLatts[2]))

get_TBL(l::Lattices.Lattice) = Lattices.NearbyUCs(l, 1)

leads = map(zip(leadlabels, LeadLatts)) do (label, latt)

	intra, inter = TBmodel.BlochHamilt_ParallPerp(latt, get_TBL, lead_hopp)

	function gf(E::Real; k=[])::Matrix{ComplexF64}

		GreensFunctions.GF_SanchoRubio(E + delta, intra(k), only(values(inter)); target="+")["+"]
	
	end


	return LayeredSystem.PrepareLead(label, latt, lead_hopp_matr, lead_hopp_matr, gf)

end 

println()


#===========================================================================#
#
# Lead lattice, spectrum and DOS 
#
#---------------------------------------------------------------------------#

for l in LeadLatts 

	k  = only(keys(l.Atoms))

	v = Lattices.Atoms_ManyUCs(l,StopStart=(0,0))

	plot_atoms(3, v; label=k, alpha=0.3, s=100)

end 

PyPlot.figure(6) 

lead_Energies = LinRange(-4.2,4.2, 200) 

lead_gf_labels = ["+","bulk"] 

DOS = zeros(2, 2,length(lead_Energies))

for (il,(l1,L1)) in enumerate(zip(leads,LeadLatts))

	H = TBmodel.Bloch_Hamilt(Lattices.NearbyUCs(L1,1);
													 argH="k",
													 lead_hopp...)

	kPoints, kTicks = Utils.PathConnect(Lattices.BrillouinZone(L1), 150; dim=2)

	spectrum = BandStructure.Diagonalize(H, kPoints; dim=2)

	PyPlot.scatter(spectrum["kLabels"], spectrum["Energy"], alpha=0.5,
								 label= only(Lattices.sublatt_labels(L1)),s=10)


	PyPlot.xlabel(spectrum["kLabel"])
	PyPlot.ylabel("Energy")

	intra, inter = TBmodel.BlochHamilt_ParallPerp(L1, get_TBL, lead_hopp)

	intra_,inter_ = intra(), only(values(inter))


	il==1 && @show intra_ inter_ 

	for e0 in LinearAlgebra.eigvals(Matrix(intra()))

		PyPlot.plot([0,1.7],[e0,e0],linestyle=":",c="gray")

	end 


	for (ie,E) in enumerate(lead_Energies)
	
		gfs = GreensFunctions.GF_SanchoRubio(E + delta, intra_, inter_; 
																				 target=join(lead_gf_labels))

#		ie==1 && 
		@assert isapprox(gfs["+"],l1[:GF](E)[1])


		for (igf,gflab) in enumerate(lead_gf_labels)

			DOS[il,igf,ie] = ObservablesFromGF.DOS(gfs[gflab])
	
		end 
	
	end 

end 

DOS = Utils.Rescale(DOS, [1,1.7], [0, maximum(DOS)])  


for (il,l1) in enumerate(leads), (igf,gflab) in enumerate(lead_gf_labels)

	PyPlot.plot(DOS[il, igf, :], lead_Energies; 
							label=string(l1[:label]," ",gflab), alpha=0.5)

end 

PyPlot.xlim(0,1.7)
PyPlot.legend() 

error()

#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

														



NG = LayeredSystem.NewGeometry(rand(0,0), 
															 default_LayerAtomRels; 
															 LeadContacts=LeadContacts,
															 Leads=leads, dim=2, 
															 dev_hopp...)


@assert length(NG)==4 

LayerAtom, Slicer, LeadRels, VirtLeads = NG 



G = GreensFunctions.GF_Decimation(dev_hopp, VirtLeads, Slicer;
																		LayerAtom...)



ENERGIES = LinRange(0.001,.5,30)


Y = zeros(4, length(ENERGIES)) 

################# 

for (iE,Energy) in enumerate(ENERGIES)

	g = G(Energy)

	get_SE = GreensFunctions.SelfEn_fromGDecim(g, VirtLeads, Slicer)

	for (ic,c) in enumerate((ObservablesFromGF.CaroliConductance,
													 ObservablesFromGF.CaroliConductance2))

		cc = c(g, get_SE, "Left", "Right", 2)

		@assert imag(cc) < 1e-10 

		Y[ic,iE] = real(cc)

	end  


	ff = ObservablesFromGF.FanoFactor2(g, get_SE, "Left", "Right", 2) 


	@assert imag(ff) <1e-10 

	Y[3,iE] = real(ff)




	print("\r",round(100iE/length(ENERGIES),digits=1),"%     ")


end 

println("\r                     ")

PyPlot.figure(4)



for (iy,(y,color)) in enumerate(zip(eachrow(Y),["blue","red","green","yellow"]))

	PyPlot.plot(ENERGIES,y,label=iy,alpha=0.6,c=color)
	PyPlot.plot(-ENERGIES,y,alpha=0.6,c=color)

end 

PyPlot.legend()






