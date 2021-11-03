import myLibs: GreensFunctions, Lattices, H_Superconductor,TBmodel
import PyPlot, LinearAlgebra

PyPlot.close.(1:10)

# tests as in
# J Comput Electron (2013) 12:203–231 
# The recursive Green’s function method for graphene
# Caio H. Lewenkopf · Eduardo R. Mucciolo



# Armchair # edges, M = 360 and N = 70 (aspect ratio W/L = 5.2) 
 
NrLayers = 10#70

NrAtomsPerLayer = 5#36#0 



latt = Lattices.HoneycombLattice() 


PyPlot.figure(2) 

for n in 0:3, m in 0:3

	PyPlot.scatter(eachrow(Lattices.Atoms_ManyUCs(latt; Ns=[n,m]))...)

end 

for arrow in Lattices.eachvec(latt.LattVec) 

	PyPlot.arrow(0,0,arrow...)
	
end  
PyPlot.gca().set_aspect(1)






Lattices.Superlattice!(latt, [[1 1]; [-1 1]]) 

# latt.Atoms["A"][:,1] .+= latt.LattVec[:,2]








#layer = unit_cell*nr_subs + sublattice 




PyPlot.figure(1) 

atoms(n,m) = Lattices.Atoms_ManyUCs(latt; Ns=[n,m])


for arrow in Lattices.eachvec(latt.LattVec) 

	PyPlot.arrow(0,0,arrow...)
	
end 


#PyPlot.scatter(eachrow(atoms(0,0))...,c="k")

for n in 1:4

#	continue 

	PyPlot.scatter(eachrow(atoms(n,0))...)
	PyPlot.scatter(eachrow(atoms(0,n))...)
	
end 

for n in 0:2,m in 0:2
	
#	PyPlot.scatter(eachrow(atoms(m,n))...)
	
end 



PyPlot.gca().set_aspect(1) 

PyPlot.close.(1:2)


PyPlot.figure(3)
PyPlot.gca().set_aspect(1) 


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

	PyPlot.scatter(eachrow(AtomsOfLayer(layer))...,
#								 label=layer)
								 c=layer in 1:NrLayers ? "k" : "red")

end 








leadlabels = ["left","right"]



LeadLatts = [Lattices.Lattice(-Lattices.LattVec(L),
															"Left"=>mapreduce(AtomsOfLayer, hcat, 0:-1:-3)),
						 Lattices.Lattice(Lattices.LattVec(L),
															"Right"=>mapreduce(AtomsOfLayer, hcat, NrLayers+1:NrLayers+4))]


for l in LeadLatts 

	for (k,v) in l.Atoms 

		PyPlot.scatter(eachrow(v)...; label=k,alpha=0.3,s=100)

	end 

end 


		

PyPlot.legend()   



#leads = map(zip(leadlabels, LeadLatts, coupling_Hparam, lead_Hparam, P3)) do (
#											label, latt, coupling_HP, lead_HP, lp)
#
#	hopps = hopp_c, hopp_l = map([coupling_HP, lead_HP]) do Params 
#
#		hopp = H_Superconductor.SC_Domain(Params,	[1])
#
#		@assert hopp[:Hopping]([20.5, -9.5], [20.5, -8.5])[1]==Params.Hopping
#
#		return  hopp 
#
#	end 
#
#
#	if !isnothing(spectra) && !haskey(spectra, "Lead")
#		
#		@show Lattices.BrillouinZone(latt)
#
#		kPoints, = Utils.PathConnect(Lattices.BrillouinZone(latt), 100, dim=2)
#
#		LeadH = TBmodel.Bloch_Hamilt(Lattices.NearbyUCs(latt,1); argH="k", hopp_l...)
#
#		LeadDiag = BandStructure.Diagonalize(LeadH, kPoints, dim=2)
#
#		spectra["Lead"] = (LeadDiag["kLabels"],LeadDiag["Energy"])
#
#	end 
#	
#
#	hc, hl = map(hopps) do hopp
#
#		function hm(args...) 
#			
#			TBmodel.HoppingMatrix(args...; hopp...)
#
#		end 
#
#	end 
#	
#	
#	@assert hopp_l[:Nr_UCs]<=1 
#
#  get_TBL(L::Lattices.Lattice) = Lattices.NearbyUCs(L, hopp_l[:Nr_UCs])
#
#	intra, inter = TBmodel.BlochHamilt_ParallPerp(latt, get_TBL, hopp_l)
#
#	# intra might be a function of k,
#	# inter might be more than one matrix (m ∈ {1,2,...})
#
#	function gf(E::Number; k=[])::Matrix{ComplexF64}
#
#		GreensFunctions.GF_SanchoRubio(E, intra(k), only(values(inter));
#																	 target="+")["+"]
#		
#	end
#
#	return LayeredSystem.PrepareLead( label, latt, hc, hl, gf)
#
#end 
#
#




















