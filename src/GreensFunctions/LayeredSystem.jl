module LayeredSystem
#############################################################################
#
#		Contains the methods needed to implement the algorithm presented in 
#
#					Lima, Dusko, Lewenkopf - 	PRB 97, 165405 (2018)
#
#			"Efficient method for computing the electronic transport properties 
#								of a multiterminal system"
#
#		The two-terminal counterpart is presented in
#
#				Lewenkopf, Mucciolo - J Comput Electron 12, 203–231 (2013)
#
#						"The recursive Green’s function method for graphene"
#
#############################################################################

import ..LA 


import ..Utils, ..Graph, ..TBmodel, ..Lattices, ..ArrayOps, ..Algebra




#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#



function LeadAtomOrder(nr_at::Int=0; dim::Int,
											 LeftLead=nothing, RightLead=nothing, 
											 kwargs...)::Dict

	LNames,LSizes = begin

		local Ls = [LeftLead,RightLead]

		local I = findall(!isnothing,Ls)

		["LeftLead","RightLead"][I],[size.(Ls[i][:head],dim) for i in I]

	end


	out = Dict()
	
	cumLSizes = ArrayOps.recursive_cumsum(LSizes, nr_at)


	for (LN,LS,cLS) in zip(LNames, LSizes, cumLSizes)
																			 
		for (uc,(n,cls)) in enumerate(zip(LS, cLS))

			out[(LN,uc)] = cls-n+1:cls

		end
		
	end

	return out

end
 



function AllAtoms(latt::Lattices.Lattice; kwargs...)::Matrix{Float64}

	AllAtoms(Lattices.PosAtoms(latt); kwargs...)

end 


function AllAtoms(Atoms::AbstractMatrix{<:Number};
									LeftLead=nothing, RightLead=nothing,
									kwargs...
									)::Matrix{Float64}

	Leads = filter(!isnothing, (LeftLead, RightLead)) 

	isempty(Leads) && return Atoms

	C(list) = cat(list...; dims=kwargs[:dim]) 

	return C((Atoms, C.(getindex.(Leads,:head))...))

end




#===========================================================================#
#
#	Groups the atoms into several layers, such that
#					- all atoms connected to leads are in layer 1
#					- for n>1, layer n contains all neighbors of atoms in layer (n-1)
#
#---------------------------------------------------------------------------#

function Distribute_Atoms(Atoms::AbstractMatrix{<:Number}, 
													isbond::Function, 
													LeadContacts::AbstractVector{<:AbstractVector{Int}}; 
													dim::Int)::Dict{Symbol, Any}
			"""	
	Atoms::Array{Float64,2} -> positions of atoms in the scattering region

isbond(a,b)::Function -> says whether there's a bond between a and b
									(possibly much faster than computing the Hamiltonian)

LeadContacts

		"""



		# -------- find layers iteratively ----------- # 


	function find_layers_iter(layer::Int, 
								 candidates::AbstractVector{Int}, 
								 out::AbstractDict{Int,Int}=Dict{Int64,Int64}()
								 )::Tuple{Int,Dict{Int,Int}}
	
		if isempty(candidates) || size(Atoms,dim)==out.count
		
			Set(keys(out))==Set(axes(Atoms,dim)) && return layer-1,out
	
			error("There are unallocated atoms.")
	
		end
	
		merge!(out, Dict(c=>layer for c in candidates))
	
		candidates = filter(setdiff(axes(Atoms,dim), keys(out))) do i
	
							any(candidates) do j 
	
								isbond(selectdim(Atoms, dim, j), selectdim(Atoms, dim, i))
	
							end 
		end 
	
	
		return find_layers_iter(layer+1, candidates, out)
	
	end
	
	
	N,D = find_layers_iter(1, vcat(LeadContacts...) |> lc -> isempty(lc) ? [1] : lc)
	
	
	# -------- sanity check and return -------------------- #
	
	
	LayerOfAtom, IndsAtomsOfLayer = Utils.FindPartners(D, sortfirst=true)
	
	
	
	
	return Dict(
	
			:NrLayers=> N,
	
			:LayerOfAtom => LayerOfAtom,
	
			:IndsAtomsOfLayer => IndsAtomsOfLayer,
	
			:AtomsOfLayer => L->selectdim(Atoms, dim, IndsAtomsOfLayer(L)),
	
						)
	
	
end


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#




function get_LeadContacts(Atoms; Leads=[], isBond=nothing,
																			LeadContacts=nothing, 
																			kwargs...)::Vector{Vector{Int}}


	isempty(Leads) && return [Int[]] 


	if LeadContacts isa AbstractVector  

		if length(Leads)==length(LeadContacts) 

			if all(!isempty(lc) && Utils.isList(lc,Int) for lc in LeadContacts)

				return LeadContacts 
			
			end  

		end 
		
	end  


	return [get_LeadContacts(Atoms, L[:head][1], isBond; kwargs...) for L in Leads]

end


function get_LeadContacts(device, lead,
													bond_len::Real; dim::Int, kwargs...
												 )::Vector{Int}

	get_LeadContacts(device, lead,
									 EuclDistEquals(bond_len; dim=dim, kwargs...))

end


function get_LeadContacts(device,
													lead::Lattices.Lattice,
													isBond::Function;
												 kwargs...)::Vector{Int} 

	get_LeadContacts(device, Lattices.PosAtoms(lead), isBond)

end 


function get_LeadContacts(device::AbstractDict,
													lead::AbstractMatrix,
													isBond::Function;
												 kwargs...)::Vector{Int}

	Utils.flatmap(1:device[:NrLayers]) do layer 

		i_relative = get_LeadContacts(device[:AtomsOfLayer](layer), lead, isBond)

		return IndsAtomsOfLayer(layer)[i_relative]

	end 

end 





function get_LeadContacts(device::Union{Lattices.Lattice, AbstractMatrix},
													lead::AbstractMatrix,
													isBond::Function;
												 kwargs...)::Vector{Int} 

	findall(any.(eachcol(isBond(lead, Lattices.PosAtoms(device)))))
														 
end 



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#



function LayerAtomRels_(Atoms::AbstractMatrix{<:Real}, 
												(LayerAtom,LeadContacts)::Tuple{AbstractDict,Any};
												kwargs...
												)::Tuple{Dict{Symbol,Any}, Vector{Vector{Int}}}

	@assert !haskey(kwargs, :LeadContacts) "Duplicate arg/kwarg"

	LayerAtomRels_(Atoms, LayerAtom; kwargs..., LeadContacts=LeadContacts)

end


function LayerAtomRels_(Atoms::AbstractMatrix{<:Real}, 
												LayerAtom::AbstractDict;
											 kwargs...
												)::Tuple{Dict{Symbol,Any}, Vector{Vector{Int}}}

	LeadContacts = get_LeadContacts(Atoms; kwargs...)

	if Check_AtomToLayer(LeadContacts; LayerAtom...)
	
		return LayerAtom, LeadContacts

	else 

		return LayerAtomRels_(Atoms, "forced";
				 										LeadContacts=LeadContacts, kwargs...)
	end

end


function LayerAtomRels_(Atoms, LayerAtom::AbstractString;
												kwargs...
												)::Tuple{Dict{Symbol,Any}, Vector{Vector{Int}}}

	@assert !haskey(kwargs,:get_leadcontacts) "Obsolete kwarg" 

	LayerAtomRels_(Atoms, Val(Symbol(lowercase(LayerAtom))); kwargs...)

end 

																#	all atoms belong to the same layer 
																
function LayerAtomRels_(Atoms::AbstractMatrix, ::Val{:trivial}; 
												kwargs...
												)::Tuple{Dict{Symbol,Any}, Vector{Vector{Int}}}

	LayerAtomRels_(([Atoms],zeros(1,1),1), Val(:sublatt); kwargs...) 

end 


function LayerAtomRels_(Atoms::AbstractMatrix, ::Val{:forced}; 
												dim::Int, isBond::Function,
												kwargs...
												)::Tuple{Dict{Symbol,Any}, Vector{Vector{Int}}}

	LeadContacts = get_LeadContacts(Atoms; kwargs...)

	out = Distribute_Atoms(Atoms, isBond, LeadContacts; dim=dim)

	return out, LeadContacts

end



function LayerAtomRels_((latt,nr_layers)::Tuple{Lattices.Lattice, Int},
												method::Val{:sublatt};
											kwargs...
											)::Tuple{Dict{Symbol,Any}, Vector{Vector{Int}}}

	LayerAtomRels_((collect(values(latt.Atoms)), 
									Lattices.LattVec(latt), nr_layers),
								 method;
								 dim=Lattices.VECTOR_STORE_DIM, kwargs...)

end 



function LayerAtomRels_((at_uc, R_, nr_layers)::Tuple{AbstractVector{<:AbstractMatrix}, <:AbstractVecOrMat, Int},
												::Val{:sublatt}; 
												dim::Int, test::Bool=false,
												kwargs...
												)::Tuple{Dict{Symbol,Any}, Vector{Vector{Int}}}

	R = Utils.VecAsMat(R_, dim) 

	@assert size(R,dim)==1

	n = length(at_uc)


	uc(layer::Int)::Int = (n-1 + layer - div(layer,n)*n)%n +1 
 
	nr_at_uc = size.(at_uc, dim)

	nr_at(layer::Int) = nr_at_uc[uc(layer)]


	function AtomsOfLayer(layer::Int)::Matrix{Float64}

#		@assert 1<=layer<=nr_layers 

		return R * floor((layer-1)/n) .+ at_uc[uc(layer)]

	end 

	
	function IndsAtomsOfLayer(layer::Int)::Vector{Int}

		@assert 1<=layer<=nr_layers 

		s = 0

		for l in 1:layer-1 

			s += nr_at(l)

		end 

		return UnitRange(s+1, s+ nr_at(layer))

	end 
	
	
	function LayerOfAtom(atom::Int)::Int 

		s = 0 

		for layer in 1:nr_layers#Int(1e10)

			na = nr_at(layer)


			1<=atom-s<=na && return layer 

			s+=na

		end 

		error("Couldn't find atom in $nr_layers layers")
	
	end 





	if test 

		for l in 1:nr_layers
			
			I = IndsAtomsOfLayer(l)
		
			@assert size(AtomsOfLayer(l),dim) == length(I)
	
			for a in I 
		
				@assert LayerOfAtom(a)==l
		
			end 

		end

	end 
	
	out = Dict{Symbol,Any}(

			:NrLayers=> nr_layers,

			:LayerOfAtom => LayerOfAtom,

			:IndsAtomsOfLayer => IndsAtomsOfLayer,

			:AtomsOfLayer => AtomsOfLayer

			)


	return (out, get_LeadContacts(out; kwargs...))


end 




#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#



function LayerAtomRels(latt::Lattices.Lattice, LayerAtom_; kwargs...
											 )::Tuple{Dict{Symbol,Any}, Vector{Vector{Int}}}

	LayerAtomRels(Lattices.PosAtoms(latt), LayerAtom_; kwargs...)

end 


function LayerAtomRels(Atoms::AbstractMatrix, LayerAtom_; 
											 kwargs...
												)::Tuple{Dict{Symbol,Any}, Vector{Vector{Int}}}

	@assert !haskey(kwargs,:get_leadcontacts) "Obsolete kwarg" 

	(LayerAtom, LeadContacts) = LayerAtomRels_(Atoms, LayerAtom_; kwargs...)

	PlotLayerAtoms_asGraph(Atoms, LayerAtom; 
												 kwargs..., LeadContacts=LeadContacts)

	return (LayerAtom, LeadContacts)

#	return get_leadcontacts ? out : LayerAtom

end


function LayerAtomRels(Atoms::Tuple, args...; kwargs...
											 )::Tuple{Dict{Symbol,Any}, Vector{Vector{Int}}}

	LayerAtomRels_(Atoms, args...; kwargs...)

end

#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function PlotLayerAtoms_asGraph(Atoms, LayerAtom;
																isBond=nothing, dim::Int,
																Leads=[], LeadContacts=nothing,
																graph_fname="",
																kwargs...)::Nothing

	isempty(graph_fname) | isnothing(Atoms) && return 

	isBond::Function 

	LeadContacts = get_LeadContacts(Atoms; Leads=Leads, isBond=isBond,
																	LeadContacts=LeadsContacts,dim=dim)

	colorrule(i) = 1<=i<=size(Atoms,dim) ? LayerAtom[:LayerOfAtom](i) : 0

	l_atoms = [cat(L[:head]...,dims=dim) for L in Leads]

	a_label(i) = string(i, i in vcat(LeadContacts...) ? "*" : "") 

	l_labels = [repeat([L[:label]],s) for (L,s) in zip(Leads,size.(l_atoms,dim))]

	labels = map(string, vcat(a_label.(axes(Atoms,dim)), l_labels...))

	Graph.PlotAtoms_asGraph(cat(Atoms,l_atoms..., dims=dim), isBond;
														colorrule = colorrule, 
														nodelabel = i->labels[i],
														fname = graph_fname)

	return 

end






#===========================================================================#
#
#	Leads and Lattice together, possibly only lattice. 
#
#---------------------------------------------------------------------------#


function NewGeometry(args...; Leads=[], nr_orb=nothing, kwargs...)

	LayerAtom, LeadContacts = LayerAtomRels(args...; Leads=Leads, kwargs...)


	VirtLeads, LeadRels = Distribute_Leads(Leads, LeadContacts; nr_orb=nr_orb, LayerAtom..., kwargs...)
	
	if isnothing(nr_orb) || (
													!isempty(Leads) && !haskey(Leads[1],:intracell)
																																			)
		return LayerAtom, LeadRels, VirtLeads 

	end
	
	Slicer = LeadLayerSlicer(;LeadRels..., LayerAtom..., VirtLeads...,
																				 				nr_orb=nr_orb)
	
#	return LayerAtom, delete!(LeadRels,:LeadSlicer), VirtLeads, Slicer 

	return LayerAtom, Slicer, delete!(LeadRels,:LeadSlicer), VirtLeads 




end







#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#




function LayerSlicer(;LayerOfAtom::Function, 
										 IndsAtomsOfLayer::Function, 
										 nr_orb::Int, kwargs...)

	T = Tuple{Tuple{AbstractString,Int},
						Tuple{Union{Colon, Vector{Int}}}}

	function slicer(name::Symbol, index::Int64)::Union{T,Nothing}

		slicer(string(name), index)

	end 


	function slicer(name::String, index::Int64)::Union{T, Nothing}


		name == "Layer" && return (name,index),(Colon(),)

		name == "Atom" || return nothing #error("Type '$name' not understood")

		layer = LayerOfAtom(index)

		isnothing(layer) && return nothing 

		atoms = IndsAtomsOfLayer(layer)

		isnothing(atoms) && return nothing

		length(atoms) == 1 && return ("Layer",layer),(Colon(),)
		
		return ("Layer",layer), (TBmodel.Hamilt_indices(collect(1:nr_orb),
															indexin(index,atoms)[1], nr_orb)
														,)
	end


	function slicer(name::Union{String,Symbol}, index::Int64, args...)#::Union{Nothing,Tuple{Tuple{String,Int,String,Int}, Tuple{Union{Colon, Vector{Int}}, Union{Colon, Vector{Int}}} 


		out1 = slicer(name, index)

		isnothing(out1) && return nothing 

		@assert out1 isa T 
		
		out2 = slicer(args...)


		isnothing(out2) && return nothing 

		@assert out2 isa T 

		(ni1,slice1),(ni2,slice2) = out1,out2

		

		return (ni1...,ni2...), (slice1...,slice2...)
	
	end 

	return slicer
		
end

#===========================================================================#
#
#	Sanity check for the atom <--> layer relationship:
#			any lead should be couple to a single & terminal layer
#
#---------------------------------------------------------------------------#


function Check_AtomToLayer(LeadContacts=[]; kwargs...)::Bool
	
	for key in [:NrLayers,:LayerOfAtom,:IndsAtomsOfLayer,:AtomsOfLayer]
		
		haskey(kwargs, key) || return false

	end

	N,LayerOfAtom = kwargs[:NrLayers], kwargs[:LayerOfAtom]

	for LC in LeadContacts

		contacts = LayerOfAtom.(LC)

		!any(boundary -> all(contacts.==boundary),[1,N]) && return false

	end

	return true

end




#==========================================================================#
#
#	Combine several (disjoint) leads into a single one 
#
#---------------------------------------------------------------------------#

function Combine_Leads(leads, atoms::AbstractMatrix{Float64}, label; 
											  dim::Int, kwargs...)

	isempty(leads) && return nothing

	coupling(l) = l[:coupling](l[:head][1],atoms)

	hamilt = haskey(leads[1],:coupling) & haskey(leads[1],:intracell)

	if length(leads)==1

		lattice_out = Dict( :label =>	label, :head => leads[1][:head])

		!hamilt && return lattice_out,nothing

		return merge(lattice_out, Dict(

			:coupling => coupling(leads[1]),

			(k=>leads[1][k] for k in [:intracell,:intercell,:GF])...
			
			)), size.(leads[1][:intracell],dim)

	end


					#	helper function
			
	f(k) = [l[k] for l in leads]

	f(k,j) = [item[min(j,end)] for item in f(k)]

#	f(k,j,E) = [item(E) for item in f(k,j)]

	nr_ucs = maximum(length, f(:head))


	lattice_out = Dict(

		:label =>	label,

		:head => map(1:nr_ucs) do j cat(f(:head,j)...,dims=dim) end,
		)

	!hamilt && return lattice_out,nothing


	NewLead = merge(lattice_out, Dict(

		:coupling => vcat(coupling.(leads)...),
																		
		:intracell => map(1:nr_ucs) do j ArrayOps.BlkDiag(f(:intracell,j)) end,

		:intercell => map(1:nr_ucs) do j ArrayOps.BlkDiag(f(:intercell,j)) end,

		:GF => function total_GF(E::Number)::Vector{Matrix{ComplexF64}}

							gfs_at_E = [l[:GF](E) for l in leads]	# one single evaluation!
							
							return map(1:nr_ucs) do j 		

								ArrayOps.BlkDiag(q[min(j,end)] for q in gfs_at_E)

							end

						end
		
						))

	subsizes = map(1:nr_ucs) do j size.(f(:intracell,j),dim) end

	return NewLead,subsizes

end




#return VirtLeads,LayerAtom,LeadRels



#===========================================================================#
#
# Groups the user-given leads into two "virtual" leads: LeftLead & RightLead
#
#---------------------------------------------------------------------------#


function Distribute_Leads(
							Leads, 
							LeadContacts::AbstractVector;#{<:AbstractVector{<:Int}}; 
							NrLayers::Int, LayerOfAtom::Function, 
							AtomsOfLayer::Function, 
							IndsAtomsOfLayer::Function, 
							nr_orb=nothing, kwargs...)

	isempty(Leads) && return Dict(),Dict{Symbol,Function}()


	lead_distrib = Dict(("LeftLead",1)=>[])

	if NrLayers!=1

		lead_distrib[("RightLead",NrLayers)] = []

	end


	
	for (iL,LC) in enumerate(LeadContacts)

		isempty(LC) && error("Lead does not couple to any atoms")
		
		layers = LayerOfAtom.(LC)

		for (side,n) in keys(lead_distrib)

			all(isequal(n), layers) && push!(lead_distrib[(side,n)], iL)

		end	
	end

	if sum(length, values(lead_distrib)) != length(Leads)

		error("A lead is not coupled to terminal layers.")

	end




	VirtLeads, LeadSizes = Dict(), Dict()

	for ((side,n),i) in filter!(p->!isempty(p.second),lead_distrib)

		out = Combine_Leads(Leads[i],	AtomsOfLayer(n), side; kwargs...)

		VirtLeads[Symbol(side)] = out[1]
		
		if !isnothing(out[2])
			LeadSizes[side] = out[2]
		end

	end

	
	SideOfLead, LeadsOfSide_ = Utils.FindPartners([string(Leads[i][:label])=>string(side) for ((side,n),I) in lead_distrib for i in I ], sortfirst=string)



	LeadsOfSide(s) = LeadsOfSide_(string(s))

	length(VirtLeads) > length(LeadSizes) && return (VirtLeads, 
								
								Dict{Symbol,Function}(	
											:SideOfLead => SideOfLead, 
											:LeadsOfSide => LeadsOfSide,
									#		:LeadSlicer => nothing,

										))


	nr_at = sum(length∘IndsAtomsOfLayer, 1:NrLayers)

	AtomsInLead,LeadOfAtom = Utils.FindPartners(LeadAtomOrder(nr_at; kwargs..., VirtLeads...))



	slicer(name::Symbol,index::Int64) = slicer(string(name),index)

	function slicer(name::String,index::Int64)


		name in ["LeftLead","RightLead"] && return (name,index),(Colon(),)

		if name=="Atom"

			lead = LeadOfAtom(index)
		
			isnothing(lead) && return nothing


			out1, (slice1, ) = slicer(lead[1]...)


			i_at = indexin(index, AtomsInLead(lead[1]))[1]

			slice2 = TBmodel.Hamilt_indices(1:nr_orb, i_at, nr_orb)


			slice1 isa Colon && return out1, (slice2, )


			return out1, (slice1[slice2],)

		end


		side = SideOfLead(name)

		allleads = LeadsOfSide(side)


		length(allleads)==1 && return (side,index),(Colon(),)

	 	lead = findfirst(isequal(name),allleads)



		R = Utils.sepLengths_cumulRanges(LeadSizes[side][min(end,index)], lead) 

		return (side,index), (R,)


	end

	
	return VirtLeads, Dict{Symbol,Function}(
													:SideOfLead => SideOfLead, 
													:LeadsOfSide => LeadsOfSide,
													:LeadSlicer => slicer
												)
end




#===========================================================================#
#
#		Brings together the (layer <--> atom) and (real lead <--> virtual lead)
#					relations into a single translator function;
#					
#			↪	 also gives the sector corresponding to the desired atom/real lead
#
#---------------------------------------------------------------------------#

function LeadLayerSlicer(;NrLayers, IndsAtomsOfLayer, LeadSlicer=nothing, nr_orb=nothing, kwargs...)

	layer_slicer = LayerSlicer(;IndsAtomsOfLayer=IndsAtomsOfLayer, nr_orb=nr_orb, kwargs...)
														
#kwargs: NrLayers, LayerOfAtom, IndsAtomsOfLayer, AtomsOfLayer, nr_orb, SideOfLead, LeadsOfSide, LeadSlicer


	function slicer(name::Union{String,Symbol}, index::Int64, args...)

		out1 = slicer(string(name), index)


		isnothing(out1) && return nothing 

		out2 = slicer(args...)



		isnothing(out2) && return nothing 

		(ni1,slice1),(ni2,slice2) = out1,out2




		return (ni1...,ni2...),(slice1...,slice2...)
	
	end





	function slicer(name_::Union{String,Symbol}, index::Int64)

		name = string(name_)

		out = layer_slicer(name, index)
		
		!isnothing(out) && return out
								# isnothing(out) means the atom is in one of the leads

		return LeadSlicer(name, index)

	end


	return slicer

end




#function LeadAtomSlicer(IndsAtomsOfLayer,NrLayers,LeadSlicer,nr_orb,VirtLeads...)



#nd 
#
#===========================================================================#
#
#	
#
#---------------------------------------------------------------------------#

function get_hoppings(dev_hopp::Function, 
											LeadLayerSlicer::Function, 
										 VirtLeads::AbstractDict)::Function

"""
	input: 
				i and j, atom indices, including leads
				Ri and Rj, the corresponding positions
	______

	output: hopping matrix between sites i and j
			
"""


# intercell is given in the u direction, uc -> uc+1 
#	coupling from lead to the system

#(master_lead, unit_cell),atom_slice 	= LeadLayerSlicer(lead_name, unit_cell)
#																= LeadLayerSlicer("Atom", atom_index)

	function get_h((i,j)::Tuple{Int,Int}, Rij=nothing)

		(K1, u1, K2, u2), slice = LeadLayerSlicer("Atom",i,"Atom",j)

		if all(isequal("Layer"),[K1,K2]) 
			
			!isnothing(Rij) && return dev_hopp(Rij...)

			error("Provide (Ri,Rj)!")

		end

		T = typeof(first(keys(VirtLeads)))
	
		K2=="Layer" && return VirtLeads[T(K1)][:coupling][slice...]
		
		K1=="Layer" && return VirtLeads[T(K2)][:coupling]'[slice...]
	

		K1==K2 || error("Lead to lead hopping?")
	

		h = VirtLeads[T(K1)][u1==u2 ? :intracell : :intercell][min(u1,u2,end)]
		
		return (u1<=u2 ? h : h')[slice...] 

	end

end


#===========================================================================#
#
#		Verify that the lead dictionary has all necessary information
#
#---------------------------------------------------------------------------#

#function CheckLead(Lead=nothing,label=nothing)
#	
#	isnothing(Lead) && return nothing
#
#	if !all(haskey.([Lead],[:intracell,:intercell,:GF,:coupling,:head]))
#		error("Incomplete lead information")
#	end
#
#	get!(Lead,:label,label)::String
#	
#	Lead[:coupling]::AbstractArray
#
#	[a::AbstractArray for k in [:intracell,:intercell] for a in Lead[k]]
#
#	typeof(Lead[:GF]) <: Function || error("Provide lead GF as a function!")
#
#  return Lead
#
#end



#===========================================================================#
#
#		Map a layered TB system to a linear graph (directed path, graph)
#				
#				Store on the graph only the unique Hamiltonians (no H.C.)
#
#---------------------------------------------------------------------------#


function LayeredSystem_toGraph(HoppMatr::Function, NrLayers::Int;
															 LeftLead=nothing, RightLead=nothing)
					"""			 
	HoppMatr(unit_cell_n,unit_cell_m) = Hamiltonian between the two layers

	NrLayers -> upper bound of n,m  (lower bound is 1)

	If any leads are specified, they must contain
		-	an identifier :label => label::String
		-	coupling matrix :coupling => U::AbstractArray
		-	hopping matrix :intercell => H::AbstractArray
		-	the GF :GF => g::Function

				"""


#	LeftLead  = CheckLead(LeftLead,"LeftLead")
#	RightLead = CheckLead(RightLead,"RightLead")

	Leads = filter(!isnothing, [LeftLead,RightLead])


	NrLeadUCs = if isempty(Leads) 0 else

									maximum([length(L[:intracell]) for L in Leads])+1

							end

#	NrLeadUCs = 5

  g = Graph.MetaDiPath(NrLayers)

	Graph.set_props!(g,Dict(:NrLayers=>NrLayers,
													:UCsLeads=>NrLeadUCs,
													:LeadLabels=>[Lead[:label] for Lead in Leads]
											 ))

	for i in 1:NrLayers
		Graph.set_props!(g,i,Dict(:type=>"Layer",
														 	:name=>("Layer",i),
														 	:H=>HoppMatr(i))
														)
	
		i>1 && Graph.set_prop!(g, i-1, i, :H, HoppMatr(i-1,i))


	end


	for Lead in Leads,i in 1:NrLeadUCs

			Graph.add_vertex!(g,Dict(
							:type=>"VirtualLead",
							:name=>(Lead[:label],i),
							:H => Lead[:intracell][min(i,end)],
							(i>1 ? [] : [:GFf => Lead[:GF],])... 
					# the GF for all lead unit cells will be stored in the lead head. 
												 ))
	end

  node = Graph.node_by_prop!(g,:name)


	
	if !isnothing(RightLead)

		Graph.add_edge!(g,NrLayers,
											node(RightLead[:label],1),
											:H, RightLead[:coupling]')

		for i in 1:NrLeadUCs-1

				Graph.add_edge!(g,node(RightLead[:label],i),
													node(RightLead[:label],i+1),
													:H, RightLead[:intercell][min(i,end)]
											)
			end # intercell is given in the u direction, as here: i -> i+1 
	end

	if !isnothing(LeftLead)

		Graph.add_edge!(g,node(LeftLead[:label],1),
											1,
											:H, LeftLead[:coupling])

		for i in 1:NrLeadUCs-1

				Graph.add_edge!(g,node(LeftLead[:label],i+1),
													node(LeftLead[:label],i),
													:H, LeftLead[:intercell][min(i,end)]'
											)
		end # intercell is given in the u direction: opposite to i+1->i 
	end


	return g

end




#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function Plot_Graph(plot_graph::Nothing=nothing, args...)

end 

function Plot_Graph(fname::AbstractString, args...)

	Plot_Graph((fname, nothing), args...)

end 

function Plot_Graph((fname,IndsAtomsOfLayer)::Tuple, NrLayers::Int, g)

	indomain(i) = 1<=i<=NrLayers

	function nodelabel(i::Int)::String
		
		n = join(g[i,:name]," ")
		
		isnothing(IndsAtomsOfLayer) && return n 

		indomain(i) || return n 
	
		return n*" ("*join(Utils.IdentifyRanges(IndsAtomsOfLayer(i))," ")*")"

	end 


	Graph.Plot_Graph(Graph.SimpleDiGraph(g),
									 fname=fname,
									 nodelabel=nodelabel,
									 colorrule=indomain)
	
end 



function LayeredSystem_toGraph(NrLayers::Int;
															 LeftLead=nothing, RightLead=nothing)
					"""			 
	NrLayers -> upper bound of n,m  (lower bound is 1)

	If any leads are specified, they must contain
		-	an identifier :label => label::String

				"""


#	LeftLead  = CheckLead(LeftLead,"LeftLead")
#	RightLead = CheckLead(RightLead,"RightLead")

	Leads = filter(!isnothing, [LeftLead,RightLead])

	NrLeadUCs = 1

  g = Graph.MetaDiPath(NrLayers)

	Graph.set_props!(g,Dict(:NrLayers=>NrLayers,
													:UCsLeads=>NrLeadUCs,
													:LeadLabels=>[Lead[:label] for Lead in Leads]
											 ))

	for i in 1:NrLayers

		Graph.set_props!(g,i,Dict(:type=>"Layer",
														 	:name=>("Layer",i),
															)
														)
	end


	for Lead in Leads, i in 1:NrLeadUCs

			Graph.add_vertex!(g,Dict(
							:type=>"VirtualLead",
							:name=>(Lead[:label],i),
												 ))
	end

  node = Graph.node_by_prop!(g,:name)


	
	if !isnothing(RightLead)

		Graph.add_edge!(g, NrLayers,
										node(RightLead[:label],1))

		for i in 1:NrLeadUCs-1

				Graph.add_edge!(g,node(RightLead[:label], i),
													node(RightLead[:label], i+1),
											)
			end 
	end

	if !isnothing(LeftLead)

		Graph.add_edge!(g,node(LeftLead[:label], 1),
											1
											)

		for i in 1:NrLeadUCs-1

				Graph.add_edge!(g, node(LeftLead[:label],i+1),
													 node(LeftLead[:label],i),
											)
		end # intercell is given in the u direction: opposite to i+1->i 
	end


	return g

end



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function get_Bonds(len_or_bond;
									 AtomsOfLayer::Function, 
									 IndsAtomsOfLayer::Function,
									 NrLayers::Int,
									kwargs...
									)::Vector{NTuple{2,Int}}

	Algebra.get_Bonds(AtomsOfLayer, 
										len_or_bond,
										IndsAtomsOfLayer,
										[(l1,l2) for l1=1:NrLayers for l2=l1:min(l1+1,NrLayers)];
										kwargs...)

end


function bondRs_fromInds(inds_bonds::AbstractVector{NTuple{2,Int}};
												 IndsAtomsOfLayer::Function,
												 AtomsOfLayer::Function,
												 NrLayers::Int,
												 kwargs...
												 )::Vector{Vector{Vector{<:Real}}}

	Algebra.bondRs_fromInds(inds_bonds, AtomsOfLayer, IndsAtomsOfLayer,
													[(l1,l2) for l1=1:NrLayers for l2=l1:min(l1+1,NrLayers)]; 
													kwargs...)

end 

#	function atoms(k) 
#
#			L = LayerOfAtom(k)
#
#			i = only(indexin(k,IndsAtomsOfLayer(L)))
#
#			return selectdim(AtomsOfLayer(L), dim, i)
#
#	end  
#
#
#
#
#
#	f1 = Utils.sel(atoms1,dim) 
#
#	f1(i) = atom_i 
#
#
#
#	map(inds_bonds) do (a,b)
#
#		map([a,b]) do k 
#
#			L = LayerOfAtom(k)
#
#			i = only(indexin(k,IndsAtomsOfLayer(L)))
#
#			return collect(selectdim(AtomsOfLayer(L), dim, i))
#
#		end 
#
#	end 
#				
#end 



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#



function plot_layers(ax_or_n; max_nr_layers::Int=30, 
										 AtomsOfLayer::Function,
										 NrLayers::Int,
										 IndsAtomsOfLayer=nothing,
										 LayerOfAtom=nothing,
										 label=nothing,
										 kwargs...)

	for layer in 1:NrLayers
	
	
	
		kw = (isnothing(label) || layer>1) ? () : Dict(:label=>"label")

		color = isodd(layer) ? "black" : "gray" 
	
		Lattices.plot_atoms(AtomsOfLayer(layer), ax_or_n;
							 color = layer>max_nr_layers ? "red" : color,
							 Utils.dict_diff(kwargs,:c,:color)...,
							 kw...)
	
		layer>max_nr_layers && break 

	end 

end 








































#############################################################################
end
