module GreensFunctions
#############################################################################

import ..LA, ..ArrayOps
#import MetaGraphs:MetaDiGraph

import ..Utils, ..TBmodel, ..Lattices

import ..Graph, ..LayeredSystem
import ..LayeredSystem#: get_node, islayer, node_right, node_left, get_graphH, get_leadlabels, get_lead_GFf

#===========================================================================#
#
# Green's function using brute-force inversion
#
#---------------------------------------------------------------------------#



""" gives the GF at energy E of an isolated system S with Hamiltonian H """

#GF2(E::Number,H::AbstractMatrix)::Matrix = inv(Matrix(E*LA.I - H))

function GF(E::Number,H::AbstractMatrix)::Matrix{ComplexF64}

	LA.inv!(LA.lu!(LA.axpy!(-1,H,Matrix{ComplexF64}(E*LA.I,size(H)))))

end 


""" the system S is now coupled to another system S'
 									such that the self-energy is SE			"""

GF(E::Number,H::AbstractMatrix,SE::AbstractMatrix)::Matrix = GF(E,H+SE)


"""
	- the system S is now coupled to n = length(args) 
				other systems {S_1, S_2, ..., S_n} with known GFs g_n

	- args is a series of tuples. Element i cooresponds to system S_i:
				args[i] = (	H( S_i -> S ), g_i )  ---> as required by SelfEn
"""

GF(E::Number,H::AbstractMatrix,args...)::Matrix = GF(E,H+SelfEn(args...))


""" the Hamiltonian of S is null """
 
GF(E::Number,args...)::Matrix = GF(E, SelfEn(args...))


""" the GF of isolated S is already known, g = inv(E-H)
							with or without coupling, the total GF is, respectively:"""

GF(g::AbstractMatrix,args...)::Matrix = GF(0.0,-inv(g),SelfEn(args...))

GF(g::AbstractMatrix)::AbstractMatrix = g


#===========================================================================#
#
#	Self energy due to coupling to a system with known GF
#
#---------------------------------------------------------------------------#

	"""
	Self energy associated with the coupling to another system

	`A` := present system; `B`:= target system with GF `g`
	
	U = HoppingMatrix(B -> A)

	"""
SelfEn(U::AbstractMatrix,g::AbstractMatrix)::Matrix = *(U',g,U)

SelfEn((U,g)::NTuple{2,AbstractMatrix})::Matrix = SelfEn(U,g)
			# if one tuple is given
			
function SelfEn(arg1::NTuple{2,AbstractMatrix}, args...)::Matrix 
	#	if more than one tuple is give

	out = SelfEn(arg1)

	for arg in args 
		out += SelfEn(arg)
	end  

	return out 

end 


#===========================================================================#
#
# Decay width
#
#---------------------------------------------------------------------------#


function DecayWidth(SE::AbstractMatrix{<:Number})::LA.Hermitian  
	
	LA.Hermitian(1im*(SE-SE'))

end 

function DecayWidth(U::AbstractMatrix{<:Number}, g::AbstractMatrix{<:Number}
									 )::LA.Hermitian 
	
	DecayWidth(SelfEn(U,g))

end




#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

lead_dir(::Val{:LeftLead})::String = "left"
lead_dir(::Val{:RightLead})::String = "right"
lead_dir(L::AbstractString)::String = lead_dir(Val(Symbol(L)))


function inter_i_dir(VirtLeads::AbstractDict,
										 LeadLayerSlicer::Function,
										 k::AbstractString,
										 uc::Int=1
										 )::Tuple{Matrix, Int, String}

	(K,i),slice = LeadLayerSlicer(k,uc)
		
	inter = VirtLeads[Symbol(K)][:intercell][min(i,end)]

	return (Matrix(view(inter, slice..., slice...)'), i+1, lead_dir(K)) 

end 



function SelfEn_fromGDecim(G::Function, 
													 VirtLeads::AbstractDict, 
													 LeadLayerSlicer::Function
													 )::Function

	function self_en(k::AbstractString, uc::Vararg{Int})::Matrix

		U,i,dir =  inter_i_dir(VirtLeads, LeadLayerSlicer, k, uc...)

		return SelfEn(U, G(k, i; dir=dir))

	end

end 


function SelfEn_fromGDecim(VirtLeads::AbstractDict,
													 LeadLayerSlicer::Function,
													 k::AbstractString,
													 uc::Vararg{Int},
													 )::Function 

	U,i,dir = inter_i_dir(VirtLeads, LeadLayerSlicer, k, uc...)

	return function self_en(G::Function)::Matrix{ComplexF64}
		
		GreensFunctions.SelfEn(U, G(k, i; dir=dir))

	end 

end 



function DecayWidth_fromGDecim1(args...)::Function 

	GreensFunctions.DecayWidth∘SelfEn_fromGDecim(args...)

end 



#===========================================================================#
#
# GF as function of Energy (fE)
#
#---------------------------------------------------------------------------#




function GF_Decimation_fE(Hopping::AbstractDict{Symbol},
											 VirtLeads::AbstractDict=Dict{Symbol,Dict{Symbol,Any}}(),
											 LeadLayerSlicer=nothing;
											 NrLayers::Int, 
											 AtomsOfLayer::Function, 
											 IndsAtomsOfLayer=nothing,
#											 graph_fname::AbstractString="",
											 kwargs...)::Function

	HoppMatr(args...) = TBmodel.HoppingMatrix(AtomsOfLayer.(args)...;
																						Hopping...)

	# will be called as 'HoppMatr(layer_n, layer_m)' or just 'HoppMatr(layer_n)'

	return GF_Decimation_fE(HoppMatr, NrLayers, VirtLeads;
											 translate=LeadLayerSlicer,
#											 plot_graph=(graph_fname,IndsAtomsOfLayer),
											 kwargs...)
	

end




function GF_Decimation_fE(HoppMatr::Function, 
											 NrLayers::Int64, VL::AbstractDict...;
											 translate=nothing, #plot_graph=nothing,
											 kwargs...)::Function

	GF_Decimation_fE(prep_GF_fE(HoppMatr, NrLayers, VL...)..., translate)

end



function GF_Decimation_fE(data_H::AbstractDict{NTuple{2,Int}, 
																		<:AbstractMatrix{ComplexF64}},
											 VirtLeads::AbstractDict=Dict{Symbol,Dict{Symbol,Any}}(),
											 Slicer=nothing;
											 NrLayers::Int, 
											 AtomsOfLayer=nothing, 		# not used 
											 IndsAtomsOfLayer=nothing,# not used 
											 kwargs...)::Function 

	GF_Decimation_fE(data_H, NrLayers, Slicer, VirtLeads; kwargs...)
	
end



function GF_Decimation_fE(data_H::AbstractDict, 
													NrLayers::Int64, VL::AbstractDict; 
													kwargs...)::Function 

	GF_Decimation_fE(data_H, NrLayers, nothing, VL; kwargs...) 

end 




function GF_Decimation_fE(data_H::AbstractDict{NTuple{2,Int}, 
																		<:AbstractMatrix{ComplexF64}},
											 NrLayers::Int64, 
											 Slicer::Union{Function,Nothing}=nothing,
											 VL::AbstractDict...;
											 kwargs...
												)::Function 
	
	GF_Decimation_fE(prep_GF_fE(NrLayers, VL...; kwargs...)..., Slicer, data_H)

end 

function GF_Decimation_fE(g::MetaDiGraph, (T,R)::Tuple{DataType,Function}, 
													args...)::Function

	function get_G(E::Number)::Function 

		GF_Decimation(g, (T,R), E, args...)
	
	end 

end


function prep_GF_fE(args...; kwargs...)::Tuple{MetaDiGraph,
																							 Tuple{DataType,Function},
																							 }

	g = LayeredSystem.LayeredSystem_toGraph(args...)

	return (g, get_reEnergy(g; kwargs...))

end  




#===========================================================================#
#
# User methods 
#
#---------------------------------------------------------------------------#

function GF_Decimation(args...; kwargs...)::Function 

	GF_Decimation!(init_storage_GF(), args...; kwargs...)

end 

function GF_Decimation!(storage::NTuple{2,Dict},
												g::MetaDiGraph, E::Number, args...;
												kwargs...)::Function

	_GF_Decimation!(storage, g, E, get_reEnergy(g, E; kwargs...), args...)


end  


function GF_Decimation!(storage::NTuple{2,Dict},
												g::MetaDiGraph, (T,R)::Tuple{DataType,Function},  
												E::Number, 
												args...
												)::Function  

	_GF_Decimation!(storage, g, E::T, R(E), args...)

end
 



#===========================================================================#
#
# Internal methods 
#
#---------------------------------------------------------------------------#

include("GFD.jl")


function _GF_Decimation!(storage::NTuple{2,Dict},
												g::MetaDiGraph,
												E::Number, rE::Number,
												data_H::AbstractDict{NTuple{2,Int}, 
																		<:AbstractMatrix{ComplexF64}}
												)::Function 

	_GF_Decimation!(storage, g, E, rE, nothing, data_H)

end 

function _GF_Decimation!(storage::NTuple{2,Dict},
												g::MetaDiGraph,
												E::Number, rE::Number,
												Slicer::Union{Function,Nothing}=nothing, 
												data_H::AbstractDict{NTuple{2,Int}, 
																		<:AbstractMatrix{ComplexF64}}...
												)::Function 

	GFD.get_G!(storage, (g, data_H..., E, rE), Slicer)

end 




#GF_Decimation!(storage, prep_GF_fE(...)..., E, Slicer) 
















#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function init_storage_GF()

	(Dict{Tuple{Int,Int,String},Matrix{ComplexF64}}(),
	 Dict{String,Vector{Matrix{ComplexF64}}}()
	 )

end 

function clear_storage_GF!(D::T)::T where T<:Tuple{Vararg{Dict}}

	empty!.(D)

	GC.gc() 

	return D 

end 



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function get_reEnergy(::Val{true}, ::Val{true})::Tuple{DataType,Function}

	(Real, identity)

end 

function get_reEnergy(::Val{true}, ::Val{false})::Tuple{DataType,Function}

	(ComplexF64, real) 

end 

get_reEnergy(::Val{false})::Tuple{DataType,Function} = (Number,identity)
	
function get_reEnergy(g::MetaDiGraph; leads_have_imag=nothing, 
											#other kwargs not accepted 
											)::Tuple{DataType,Function}

	if !isempty(get_leadlabels(g)) && leads_have_imag isa Bool 

		return get_reEnergy(Val(true), Val(leads_have_imag))

	end 

	return get_reEnergy(Val(false))

end 

	
function get_reEnergy(g::MetaDiGraph, Energy::Number; kwargs...)::Number 

	T,F = get_reEnergy(g; kwargs...)

	Energy::T 

	return F(Energy)

end 


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#



#	if !isnothing(plot_graph)
#
#		graph_fname,IndsAtomsOfLayer = plot_graph
#
#		if !isempty(graph_fname)
#		
#			nodelabel2(i) = (1<=i<=NrLayers ? (" ("*join(Utils.IdentifyRanges(IndsAtomsOfLayer(i))," ")*")") : "")
#
#			nodelabel(i) = join(g[i,:name]," ")* (isnothing(IndsAtomsOfLayer) ? "" : nodelabel2(i))
#
#			Graph.Plot_Graph(Graph.SimpleDiGraph(g),fname=graph_fname,nodelabel=nodelabel,colorrule= i-> 1<=i<=NrLayers )
#
#		end
#	end

#	LayeredSystem.Plot_Graph(plot_graph, NrLayers, g_withH) 
































#===========================================================================#
#
#		Bulk, left-surface and right-surface GF, as presented in 
#
#				M. P. Lopez Sancho, J. M. Lopez Sancho, J. Rubio
#
#					J. Phys. F: Met. Phys. 15, 851-858 (1985) 
#
#			"Highly convergent schemes for the calculation of bulk and
#																surface Green functions"
#
#---------------------------------------------------------------------------#

function GF_Surface(latt::Lattices.Lattice,
										target_arg::AbstractString,
										args...
									)::Function

	only ∘ GF_Surface(latt, [target_arg], args...)


end

function GF_Surface(latt::Lattices.Lattice,
										target_arg::AbstractVector{<:AbstractString},
										hopping::AbstractDict,
										delta::Real,
										)::Function

	intra, inter = TBmodel.BlochHamilt_ParallPerp(latt, Lattices.NearbyUCs, hopping)
	
	@assert Lattices.LattDim(latt)==1 

	intra_, inter_ = intra(), only(values(inter))

	target=join(target_arg)

	return function gf(E::Real)::Vector{<:Matrix{ComplexF64}}
		
		gfs = GreensFunctions.GF_SanchoRubio(E+delta*im, intra_, inter_;
																				 target=target)

		return [gfs[t] for t in target_arg]

	end

end



function GF_Surface(latt::Lattices.Lattice,
										target_arg::AbstractVector{<:AbstractString},
										hopping::AbstractDict,
										)::Function

	intra, inter = TBmodel.BlochHamilt_ParallPerp(latt, Lattices.NearbyUCs, hopping)
	
	@assert Lattices.LattDim(latt)==1 

	intra_, inter_ = intra(), only(values(inter))
		
	
	return function gf(E::ComplexF64)::Vector{Matrix{ComplexF64}}

		@assert abs(imag(E))>1e-8 "Small imaginary part needed"

		gfs = GreensFunctions.GF_SanchoRubio(E, intra_, inter_;
																				 target=join(target_arg))

		return [gfs[t] for t in target_arg]

	end

end







function SanchoRubio_converged!(Errors::AbstractVector{<:Real},
																iter::Int, 
																matrices::Vararg{<:AbstractMatrix};
																tol::Float64=1e-8,
																kwargs...
																)::Bool

	for M in matrices 
		
		Errors[iter] += LA.norm(M, Inf)

	end  

	nr_last = 3 

	last_few = view(Errors, max(1,iter-nr_last+1):iter)


	if (iter==length(Errors) || isinf(Errors[iter]) || isnan(Errors[iter])
			) || (iter>=nr_last && all(>(1e10), last_few))

		@warn "The Sancho-Rubio algorithm did not converge after $iter iterations. The errors are:\n $(Errors[1:iter])"

		return true
		
	end 



	return iter>=nr_last && all(<(tol), last_few)


end

function SanchoRubio_recursion(
								Energy::Number,

								max_iter::Int, 
								alpha::AbstractMatrix{<:Number},
								epsBulk::AbstractMatrix{<:Number},

								beta::AbstractMatrix{<:Number}=alpha',
								epsSurf::AbstractMatrix{<:Number}=epsBulk,
								epsDualSurf::AbstractMatrix{<:Number}=epsBulk,

								Errors::AbstractVector{Float64}=zeros(max_iter),
								iter::Int=1;
								kwargs...)::NTuple{3,Matrix{ComplexF64}}


	if SanchoRubio_converged!(Errors, iter, alpha, beta; kwargs...)

		return (epsBulk, epsSurf, epsDualSurf)

	end 


  gBulk = GF(Energy, epsBulk)
																	# auxiliary vars
  agb = alpha*gBulk*beta 

	bga = beta*gBulk*alpha


  return SanchoRubio_recursion(
															 Energy,
															 max_iter,
															 
															 alpha*gBulk*alpha, # new alpha 
															 epsBulk + agb + bga, #new epsBulk 

															 beta*gBulk*beta,	# new beta 
															 epsSurf + agb,				#new epsSurf
															 epsDualSurf + bga,		#new epsDualSurf

															 Errors,
															 iter+1;			#next iteration 
															 kwargs...
															 )

end 
 


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function SanchoRubio_inplace(
								Energy::Number,
								max_iter::Int,  

								H_intracell::AbstractMatrix{<:Number},
								H_intercell::AbstractMatrix{<:Number},

								k::Int;

								kwargs...)::Matrix{ComplexF64}

	SanchoRubio_inplace(Energy, max_iter, H_intracell, H_intercell,
											k==1, k==2, k==3; kwargs...)[k]

end  




function SanchoRubio_inplace(
								Energy::Number,
								max_iter::Int,  

								H_intracell::AbstractMatrix{<:Number},
								H_intercell::AbstractMatrix{<:Number},

								e1::Bool, e2::Bool, e3::Bool;	

								kwargs...)::NTuple{3,AbstractMatrix{ComplexF64}}

	@assert LA.ishermitian(H_intracell)  H_intracell

	Errors = zeros(max_iter)




#	alpha = copy(Matrix(H_intercell))
#	beta = copy(alpha')
# eps = cat((H_intracell for i in 0:(e2|e3))..., dims=3)


	n = LA.checksquare(H_intercell)
	@assert n==LA.checksquare(H_intracell)


	alpha = Matrix{ComplexF64}(undef, n, n)
	beta = Matrix{ComplexF64}(undef, n, n)
	eps = Array{ComplexF64,3}(undef, n, n, 1+(e2|e3))
	
	ga = Matrix{ComplexF64}(undef, n, n)
	gb = Matrix{ComplexF64}(undef, n, n) 

	aux = Matrix{ComplexF64}(undef, n, n)


	copy!(alpha, H_intercell)
	copy!(beta, alpha') 

	for k=axes(eps,3)

		copy!(selectdim(eps,3,k), H_intracell) 

	end 



	gBulk = GF(Energy, selectdim(eps, 3, 1))


	for iter in 1:max_iter  

		SanchoRubio_converged!(Errors, iter, alpha, beta; kwargs...) && break  

		LA.mul!(ga, gBulk, alpha)  
		LA.mul!(aux, beta, ga)  

		selectdim(eps, 3, 1) .+= aux  #bga 


		LA.mul!(gb, gBulk, beta)  
		LA.mul!(aux, alpha, gb)

		for k=axes(eps,3) 
			
			selectdim(eps,3,k) .+= aux # agb 

		end  

		copy!(alpha, LA.mul!(aux, alpha, ga))

		copy!(beta, LA.mul!(aux, beta, gb))

		copy!(gBulk, GF(Energy, selectdim(eps, 3, 1)))

	end 


	return (gBulk, 
					e2 ? GF(Energy, selectdim(eps,3,2)) : gBulk,
					e3 ? GF(Energy, -(eachslice(eps,dims=3)...) + H_intracell) : gBulk
					)

end 

#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function understand_target_GF(target::AbstractString
															)::Tuple{Vector{String},BitVector}

	input_keys = [["bulk"],["+","plus","positive"],["-","minus","negative"]]

	desired_keys = [any(occursin.(k,[target])) for k in input_keys]

	any(desired_keys) || return understand_target_GF()

	return ["bulk","+","-"], desired_keys

end  



function understand_target_GF(::Nothing)::Tuple{Vector{String},BitVector}

	@warn "Target not understood. Returning all."

	return ["bulk","+","-"], [true,true,true] 

end 


function understand_target_GF(;target=nothing, kwargs...
														 )::Tuple{Vector{String},BitVector}

	understand_target_GF(target)

end 


	



function GF_SanchoRubio(Energy::Number, 
												H_intracell::AbstractMatrix{<:Number},
												H_intercell::AbstractMatrix{<:Number};
												max_iter::Int=50,
												kwargs...
												)::Dict{String,Matrix{ComplexF64}}

				# The system is assumed homogeneous
				# in general, H(cell_i,cell_j) and H(cell_i)


	# ------ understanding target ------ #
	
	output_keys, desired_keys = understand_target_GF(;kwargs...)


	# ----- performing the iterations ------ # 


	gs = SanchoRubio_inplace(Energy, max_iter, H_intracell, H_intercell,
													 desired_keys...; kwargs...)


	# ----- return the desired GF ------- #


	return Dict(k=>g for (k,g,d) in zip(output_keys, gs, desired_keys) if d)

end




function GF_SanchoRubio(Energy::Number, 
												H_intracell::AbstractMatrix{<:Number},
												H_intercell::AbstractMatrix{<:Number},
												target::AbstractString;
												max_iter::Int=50,
												kwargs...
												)::Matrix{ComplexF64}

	SanchoRubio_inplace(Energy, max_iter, H_intracell, H_intercell,
											only(findall(understand_target_GF(target)[2]));
											kwargs...)

end




function GF_SanchoRubio_rec(Energy::Number, 
												H_intracell::AbstractMatrix{<:Number},
												H_intercell::AbstractMatrix{<:Number};
												max_iter::Int=50,
												kwargs...
												)::Dict{String,Matrix{ComplexF64}}

				# The system is assumed homogeneous
				# in general, H(cell_i,cell_j) and H(cell_i)

#	@warn "Obsolete implementation"

	# ------ understanding target ------ #

	output_keys, desired_keys = understand_target_GF(;kwargs...)


	# ----- performing the iterations ------ # 
	
	@assert LA.ishermitian(H_intracell)

	hs = SanchoRubio_recursion(Energy, max_iter, 
														 Matrix(H_intercell), Matrix(H_intracell);
														 kwargs...)

	# ----- return the desired GF ------- #

	return Dict(k => GF(Energy, h) for (k,h,d) in zip(output_keys, hs, desired_keys) if d)

end







#===========================================================================#
#
# Helper function, prepares the lead dictionary for further operations
# 			:label and :coupling must be additionally provided
#
#---------------------------------------------------------------------------#


""" 
	- Lead -- aligned and attached to Atoms

	- BridgeAtoms::Matrix{Float64} is non-empty only for unregular lead-atom
				iterfaces. Contains additional atoms to ensure maximum contact

	- HoppMatr::Function(atoms1,atoms2) gives the Hamiltonian

	- LeadGF::Function(E) lead GF, computed elsewhere

"""

#PrepareLead(tup::Tuple)::Dict{Symbol,Any} = PrepareLead(tup...)

function PrepareLead(atoms::AbstractMatrix{<:Real})::Dict{Symbol,Any}

	Dict{Symbol,Any}(:head=>[atoms])

end 


function PrepareLead(label::AbstractString, 
										 LeadAtoms::AbstractMatrix{<:Real}
										 )::Dict{Symbol,Any}

	merge!(Dict{Symbol,Any}(:label => label), PrepareLead(LeadAtoms))
	
end  


function PrepareLead(label::AbstractString, Lead::Lattices.Lattice
										 )::Dict{Symbol,Any}

	PrepareLead(label, Lattices.PosAtoms(Lead))

end 


function PrepareLead(label::AbstractString, 
										 Lead::Union{Lattices.Lattice, AbstractMatrix{<:Real}},
										 BridgeAtoms::AbstractMatrix{<:Real}
										 )::Dict{Symbol,Any}

	merge!(vcat, PrepareLead(BridgeAtoms), PrepareLead(label, Lead)) 

end 

function PrepareLead(label::AbstractString,
										 lead_latt::Lattices.Lattice,
										 lead_hopp::AbstractDict,
										 delta::Real...,
										 )::Dict{Symbol,Any} 

	PrepareLead(label, lead_latt, lead_hopp, lead_hopp, delta...)

end 

function PrepareLead(label::AbstractString,
										 lead_latt::Lattices.Lattice,
										 coupl_hopp::AbstractDict,
										 lead_hopp::AbstractDict,
										 delta::Real...,
										 )::Dict{Symbol,Any}

	PrepareLead(label, 
							lead_latt,
							Utils.add_args_kwargs(TBmodel.HoppingMatrix; coupl_hopp...),
							Utils.add_args_kwargs(TBmodel.HoppingMatrix; lead_hopp...),
							GF_Surface(lead_latt, "+", lead_hopp, delta...))

end



function PrepareLead(label::AbstractString, 
										 Lead::Lattices.Lattice,
										 coupling::Function,
										 LeadHoppMatr::Function,
										 LeadGF::Function,
										 )::Dict{Symbol,Any}

	LeadAtoms0, LeadAtoms1 = [Lattices.Atoms_ManyUCs(Lead; Ns=n) for n=[0,1]]

	return merge!(PrepareLead(label, LeadAtoms0), Dict{Symbol,Any}(	

						:coupling => coupling,

						:intracell => [LeadHoppMatr(LeadAtoms0)],

						:intercell => [LeadHoppMatr(LeadAtoms0, LeadAtoms1)],

						:GF => function GF(E::Number)::Vector{<:Matrix} 
						
												[LeadGF(E)] 
											
									end

								)
							)

end 



function PrepareLead(label::AbstractString, 
										 Lead::Lattices.Lattice,
										 BridgeAtoms::AbstractMatrix,
										 coupling::Function,
										 LeadHoppMatr::Function,
										 LeadGF::Function,
										 )::Dict{Symbol,Any}

	BridgeIntra = LeadHoppMatr(BridgeAtoms) 

	BridgeToLead = LeadHoppMatr(BridgeAtoms, Lattices.PosAtoms(Lead))


	out = merge!(vcat,
							 Dict{Symbol,Any}(:intracell=>[BridgeIntra], 
																:intercell=>[BridgeToLead]),
							PrepareLead(label, Lead, coupling, LeadHoppMatr, LeadGF))
				
				
	out[:GF] = function gf(E::Number)::Vector{<:Matrix}
	
				g = LeadGF(E)

				return [GF(E, BridgeIntra, BridgeToLead', g), g]

		end

	return out

end
























































#############################################################################
end

