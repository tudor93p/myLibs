module GreensFunctions
#############################################################################

import ..LA, ..ArrayOps
import MetaGraphs:MetaDiGraph

import ..Utils, ..TBmodel, ..Lattices

import ..Graph, ..LayeredSystem
import ..LayeredSystem: get_node, islead, islayer, node_right, node_left, get_graphH, get_leadlabels

#===========================================================================#
#
# Green's function using brute-force inversion
#
#---------------------------------------------------------------------------#



""" gives the GF at energy E of an isolated system S with Hamiltonian H """

GF(E::Number,H::AbstractMatrix)::Matrix = inv(Matrix(E*LA.I - H))



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

	A := present system; B:= target system with GF g
	
	U = HoppingMatrix(B -> A)

	"""

SelfEn(U::AbstractMatrix,g::AbstractMatrix)::Matrix =  U'*g*U

SelfEn((U,g)::Tuple{AbstractMatrix,AbstractMatrix})::Matrix = SelfEn(U,g)
			# if one tuple is given
			
SelfEn(args...)::Matrix = sum(SelfEn.(args))
			#	if more than one tuple is give




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
#	Helper functions for dealing with the physical system stored in the graph
#
#---------------------------------------------------------------------------#


function GraphLayeredSystem_Utils(g)

	get_prop(p) = Graph.get_prop(g,p)

	get_prop(p,ns...) = Graph.get_prop(g,ns...,p)


	node0 = Graph.node_by_prop!(g,:name)



	function node(T::String,I::Int64)

		occursin("Lead",T) && return node0(T,min(get_prop(:UCsLeads),I))

		return node0(T,I) 
	end

	node((T,I)) = node(T,I)


	f() = ()
	
	f(n::Int64,args...) = (n,f(args...)...)
	
	f(T::String,I::Int64,args...) = (node(T,I),f(args...)...)
		

		
	
	
	


	islead_(n) = occursin("Lead",get_prop(:type,n))

#	islayer_(n) = get_prop(:type,n)=="Layer"


	function H(args...)::AbstractMatrix

		ns = f(args...)
		
		(length(ns) == 1 || ns[1]==ns[2]) && return get_prop(:H,ns[1])
	
		Graph.has_prop(g,ns...,:H) && return get_prop(:H,ns...)
	
		rns = reverse(ns)
	
		Graph.has_prop(g,rns...,:H) && return get_prop(:H,rns...)'

		return zeros(size(get_prop(:H,ns[1]),1),size(get_prop(:H,ns[2]),1))

	end
	
	
	
	function LeadGF_atEnergy(g_,Energy)::Function

		all_lead_gfs = Dict(ll=>get_prop(:GFf,f(ll,1)[1])(Energy) for ll in get_prop(:LeadLabels))
 
		return function get_gf(n::Int)::AbstractMatrix 

			@assert !Graph.has_prop(g_, n, :GF) 

			lead_name, lead_uc = Graph.get_prop(g_, n, :name)

			return all_lead_gfs[lead_name][min(lead_uc,end)]

		end  

	end
	
	
	
	function left_(n::Int)::Union{Nothing,Int}
	
	  ns = Graph.in_neighbors(g,n)
	
		return isempty(ns) ? nothing : only(ns)
	
	end
	
	function right_(n::Int)::Union{Nothing,Int}
	
	  ns = Graph.out_neighbors(g,n)
	
		return isempty(ns) ? nothing : only(ns)
	
	end
	
	
	function next(n::Int64, dir::String)::Union{Nothing,Int}
	
		dir == "right" && return right_(n)
	
		dir == "left" && return left_(n)
	
		error("dir should be left or right")
	
	
	end
	


	function meets_layer_(n,dir::String="both")::Bool
	
		isnothing(n) && return false

		(dir == "both" || islayer(g,n)) && return true
	
		return meets_layer_(next(n,dir), dir)
	
	end
	
	
	
	
	
	function bring_closer_(n,m) 
	
		isempty(Graph.FindPath(g,n,m)) &&	return false,left_(n),right_(m)
	
		return true,right_(n),left_(m)
	
	end
	
	function lead_extends_(n,dir)
	
		return islead(g,n) && (dir == "both" || !meets_layer_(n,dir))
	
	end



	return 	LeadGF_atEnergy,#setEnergy_LeadGF,
					islead_,
					meets_layer_,
					lead_extends_,
					H,
					next,
					bring_closer_,
					f


end	


#===========================================================================#
#
# 	Green's function using decimation method for at most two leads.
#
#		Implementation as explained by Lewenkopf+Mucciolo in
#
#		"The recursive Green’s function method for graphene"
#
#				J Comput Electron (2013) 12: 203–231
#
#---------------------------------------------------------------------------#

function GF_Decimation(Hopping::AbstractDict, 
											 VirtLeads::AbstractDict=Dict(),
											 LeadLayerSlicer=nothing;
											 NrLayers::Int, 
											 AtomsOfLayer::Function, 
											 IndsAtomsOfLayer=nothing,
											 graph_fname::AbstractString="",
											 kwargs...)

	HoppMatr(args...) = TBmodel.HoppingMatrix(AtomsOfLayer.(args)...;
																						Hopping...)

	# will be called as 'HoppMatr(layer_n, layer_m)' or just 'HoppMatr(layer_n)'



	return GF_Decimation(HoppMatr, NrLayers;
											 VirtLeads..., 
											 translate=LeadLayerSlicer,
											 plot_graph=(graph_fname,IndsAtomsOfLayer),
											 kwargs...)
	


end




function GF_Decimation(HoppMatr::Function, NrLayers::Int64; 
											 LeftLead=nothing, RightLead=nothing, 
											 translate=nothing, plot_graph=nothing,
											 kwargs...)

	g = LayeredSystem.LayeredSystem_toGraph(HoppMatr, NrLayers; 
																					LeftLead=LeftLead, 
																					RightLead=RightLead)


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

	LayeredSystem.Plot_Graph(plot_graph, NrLayers, g)

	return function gf(Energy::Number)
		
		GF_Decimation_fromGraph(Energy, g, translate; kwargs...)

	end 

end


function GF_Decimation_fromGraph(Energy::Number, g, translate=nothing;
#																 leads_have_imag::Bool=false,
																 kwargs...
																 )::Function


for l in get_leadlabels(g)


n = get_node(g,l,1) 

h = Graph.get_prop(g,n,:H) 
@assert LA.ishermitian(h)

end 


#	setEnergy_LeadGF,
	LeadGF_atEnergy,
	islead_,
	meets_layer_,
	lead_extends_,
	H,
	next,
	bring_closer_,
	node	= GraphLayeredSystem_Utils(g)


	reEnergy = if !isempty(get_leadlabels(g))&&haskey(kwargs, :leads_have_imag)

				kwargs[:leads_have_imag] ? Energy::Real : real(Energy::Complex)

							else 

								Energy

							end
	

	SemiInfLeadGF = LeadGF_atEnergy(g, Energy)


	# ----- the dictionary where the values are stored -------- #
	
	Gd = Dict{Tuple{Int64,Int64,String},Array{Complex{Float64},2}}()



	function G(n::Int,m::Int,dir::AbstractString)::AbstractMatrix
	
		n==m && islead_(n) && !meets_layer_(n,dir) && return SemiInfLeadGF(n)

							#	this should not be stored, already stored in the graph!

		n==m && return get!(Gd,(n,m,dir)) do 
																			Gnn(n,dir) end

		return get!(Gd,(n,m,dir)) do 
																Gnm(n,m,dir) end

	end


	# ----- methods to calculate the GF ------- #	

	function Gnn(n::Int, dir::AbstractString)::Matrix{ComplexF64}
	
		lead_extends_(n,dir) && return GF(SemiInfLeadGF(n),
																			coupling_toSystem(n,dir)...)
	
		return GF(reEnergy, H(n), coupling(n,dir)...) 
	
	end

	function Gnm(n::Int, m::Int, dir::AbstractString)::Matrix{ComplexF64}

		isnleft, n1, m1  = bring_closer_(n,m)
	
		for (d,test) in zip(["left","right"],[isnleft,!isnleft])
	
			if dir in ["both",d] 
	
				test && return G(n,m1,d)*H(m1,m)*G(m,m,dir)
	
				return G(n,n,dir)*H(n,n1)*G(n1,m,d)
	
			end
		end

	end



	# ---- coupling a layer left or right ----- #

	function coupling(src::Int64,dir::String)
	
		filter(!isempty, map(["right","left"]) do d
	
										dst = (dir in [d,"both"]) ? next(src,d) : nothing
	
										return isnothing(dst) ? () : (H(dst,src), G(dst,dst,d))
	
									end)
	end
	
	
	# ------- coupling for a lead ---------- #

	
	function coupling_toSystem(src,d) 

		if d=="both" 

			@info "Function called for d='both'"

		else 

			@warn "potentially wrong result" src d

		end 

	
		for dir in ["right","left"]
	
			meets_layer_(src,dir) && return coupling(src,dir)
	
		end
	
	end
	



	if	isnothing(translate)
	
		function out1(name1::AbstractString,index1::Int64,
								 	name2::AbstractString=name1,index2::Int64=index1;
									dir::AbstractString="both")::AbstractMatrix
	
			G(node(name1,index1,name2,index2)...,dir)
	
		end


		function out1((name1,index1)::Tuple{AbstractString,Int},
									(name2,index2)::Tuple{AbstractString,Int}=(name1,index1);
									kwargs1...)::AbstractMatrix

			out1(name1, index1, name2, index2)

		end 

		function out1(ni1ni2::Tuple{Tuple,Tuple}; kwargs1...)::AbstractMatrix

			out1(ni1ni2...)

		end

		return out1

	
	else


		function out2(name1::AbstractString, index1::Int64,
								 	name2::AbstractString=name1, index2::Int64=index1;
									dir::AbstractString="both")::AbstractMatrix

			n1i1n2i2,slice = translate(name1,index1,name2,index2) 

			@assert length(slice)==2

			return G(node(n1i1n2i2...)...,dir)[slice...]

		end 

		function out2((name1,index1)::Tuple{AbstractString,Int},
									(name2,index2)::Tuple{AbstractString,Int}=(name1,index1);
									kwargs1...)::AbstractMatrix

			out2(name1, index1, name2, index2)

		end 

		function out2(ni1ni2::Tuple{Tuple,Tuple}; kwargs1...)::AbstractMatrix

			out2(ni1ni2...)

		end

	end

end



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function eval_lead_GF(g::MetaDiGraph, Energy::Number 
											)::Function

	all_lead_gfs = Dict{String,Vector{Matrix{Complex}}}()

	for label in Graph.get_prop(g, :LeadLabels) 

		GF = Graph.get_prop(g, get_node(g, label, 1), :GFf)

		all_lead_gfs[label] = GF(Energy)
		
	end 

	return function get_gf(n::Int)::AbstractMatrix{ComplexF64}

		lead_name, lead_uc = Graph.get_prop(g, n, :name)

		return all_lead_gfs[lead_name][min(lead_uc,end)]

	end  

end




function get_reEnergy(Energy::T, ::Val{false})::T where T<:Number 

	Energy

end 

function get_reEnergy(Energy::T, ::Val{true}, ::Val{true})::T where T<:Real

	Energy 

end 


function get_reEnergy(Energy::Complex, ::Val{true}, ::Val{false})::Real 

	real(Energy)

end 


	
	
function get_reEnergy(g::MetaDiGraph, Energy::Number; kwargs...)::Number 

	
	if !isempty(get_leadlabels(g)) && haskey(kwargs, :leads_have_imag)

		return get_reEnergy(Energy, Val(true), Val(kwargs[:leads_have_imag]))

	end 


	return get_reEnergy(Energy, Val(false))

end 

function next_node(g::MetaDiGraph, 
									 n::Int64, 
									 dir::AbstractString)::Union{Nothing,Int}

	dir=="right" && return node_right(g, n)

	dir=="left" && return node_left(g, n)

	error("dir should be left or right")

end 

meets_layer(::MetaDiGraph, ::Nothing, ::AbstractString="")::Bool = false  

meets_layer(::MetaDiGraph, ::Int)::Bool = true  

function meets_layer(g::MetaDiGraph, n::Int, dir::AbstractString)::Bool

	(dir=="both" || islayer(g,n)) && return true

	return meets_layer(g, next_node(g,n,dir), dir)

end

function get_dir_layer(g::MetaDiGraph, n::Int)::String 

	for dir in ("right","left")

		meets_layer(g,n,dir) && return dir 
		
	end 
	
	error()

end 


function lead_extends(g::MetaDiGraph, n::Int, dir::AbstractString)::Bool

	islead(g,n) && (dir=="both" || !meets_layer(g,n,dir))

end



# ---- coupling node(s) left or right ----- #
	
function coupling_inds(g::MetaDiGraph, src::Int, 
														 dir::AbstractString
														 )::Vector{Tuple{Int,Int,String}}

	if dir in ("right","left")

		 dst = next_node(g, src, dir) 

		 return isnothing(dst) ? Tuple{Int,Int,String}[] : [(src,dst,dir)]

	end 

	return vcat((coupling_inds(g, src, d) for d=("right","left"))...)

end  


function bring_closer(g::MetaDiGraph, n::Int, m::Int
											)::Tuple{Bool,Int,Int}

	if isempty(Graph.FindPath(g,n,m))  #means m<n
		
		return (false, node_left(g,n), node_right(g,m))

	else 

		return (true, node_right(g,n), node_left(g,m))

	end

end  

function baz(g::MetaDiGraph, n::Int, m::Int, dir::AbstractString
						 )::NTuple{2,Tuple{Int,Int,String}}

	nisleft, n1, m1  = bring_closer(g, n, m)
	
#	for (d,test) in zip(["left","right"],[nisleft,!nisleft])
#	
#		if dir in ["both",d] 
#	
#			return test ? ((n,m1,d), (m,m,dir)) : ((n,n,dir), (n1,m,d)) 
#
#		end 
#	
#	end 

	if dir in ("both","left")

		return nisleft ? ((n,m1,"left"), (m,m,dir)) : ((n,n,dir), (n1,m,"left"))  

	elseif dir=="right"

		return nisleft ? ((n,n,dir), (n1,m,dir)) : ((n,m1,dir), (m,m,dir)) 

	end 

end 


function unpack_argsG(ni1ni2::Tuple{Tuple,Tuple}; kwargs1...
							)::Tuple{Tuple{String,Int,String,Int},String}

	unpack_argsG(ni1ni2...; kwargs1...)

end 
function unpack_argsG(n1i1::Tuple{<:AbstractString,Int},
											n2i2::Tuple{<:AbstractString,Int}=n1i1;
											kwargs1...
							)::Tuple{Tuple{String,Int,String,Int},String}

	unpack_argsG(n1i1..., n2i2...; kwargs1...)

end 


function unpack_argsG(name1::AbstractString, index1::Int64,
						 	name2::AbstractString=name1, index2::Int64=index1;
							dir::AbstractString="both"
							)::Tuple{Tuple{String,Int,String,Int},String}

	(name1,index1,name2,index2), dir 

end 



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function GF_Decimation_fromGraph(Energy::Number, 
																 g::MetaDiGraph,
																 data_H::Dict;
																 kwargs...
																 )::Function

	reEnergy = get_reEnergy(g, Energy; kwargs...)

	SemiInfLeadGF = eval_lead_GF(g, Energy)


	# ----- the dictionary where the values are stored -------- #
	
	Gd = Dict{Tuple{Int,Int,String},Matrix{ComplexF64}}() # trapped storage 


	function UGpair((src,dst,dir)::Tuple{Int,Int,String},
									)::Tuple{<:AbstractMatrix{ComplexF64},
													 <:AbstractMatrix{ComplexF64},
													 }

		(get_graphH(g,data_H,dst,src; zeros_missing=true), G(dst,dst,dir))

	end  


	# ----- function returned -------- # 

	function G(n::Int,m::Int,dir::AbstractString)::AbstractMatrix{ComplexF64}
	
		n==m && islead(g,n) && !meets_layer(g,n,dir) && return SemiInfLeadGF(n)

							#	this should not be stored, already stored in the graph!

		n==m && return get!(Gd,(n,m,dir)) do 
																			Gnn(n,dir) end

		return get!(Gd,(n,m,dir)) do 
																Gnm(n,m,dir) end

	end



	# ----- methods to calculate the GF ------- #	

	function Gnn(n::Int, dir::AbstractString)::Matrix{ComplexF64}

		if lead_extends(g,n,dir) 

			dir0 = get_dir_layer(g,n)

		if dir=="both" 

			@info "Function called for dir='both'"

		else 

			@warn "potentially wrong result" n dir  dir0

			@show coupling_inds(g, n, dir) 

			@show coupling_inds(g, n, dir0)

			error() 
		end  

#			dir=="both" || !meets_layer(g,n,dir) 

# n is lead 
#
# 1) Gnn(n, dir=="both"): GF(g_n, (U_coupling,G_system))  regardless of dir 
#
# 2) Gnn(n,dir=="right") 
#	lead semi-infinite in dir=="right". !meets_layer(g,n,"right")
#
#	G^R_{N+1,N+1} from Lewenkopf 
#
#	(30): G^R_{N,N} = (E - h_N - U_{N,N+1}g_R U_{N+1,N}
#
#	G^R_{Lead R, uc 1} = (E - intracell - inter*g*inter') # no layer etc
#
#	lci = coupling_inds(g, n, dir=get_dir_layer(g,n)="left")
#	lci = (n,n_left,"left"), n_left either lead_uc-1 or layer N
#
#	UGpair(right_lead, layer_n, "left") 
# lci=coupling_inds(g, n, dir)
#
			return GF(SemiInfLeadGF(n),
								(UGpair(lci) for lci=coupling_inds(g, n, dir0))...)

		end 
																										 
	
	
		return GF(reEnergy, 
							get_graphH(g, data_H, n; zeros_missing=true),
							(UGpair(lci) for lci=coupling_inds(g, n, dir))...
							)
	end


	function Gnm(n::Int, m::Int, dir::AbstractString)::Matrix{ComplexF64}

		(a1,b1,d1),(a2,b2,d2) = baz(g, n, m, dir) 

		return *(G(a1,b1,d1),
						 get_graphH(g, data_H, b1, a2; zeros_missing=true),
						 G(a2,b2,d2),
						)

	end 


	return function out1(args...; kwargs1...)::AbstractMatrix{ComplexF64}
	
		(n1,i1,n2,i2), dir = unpack_argsG(args...; kwargs1...) 

		return G(get_node(g,n1,i1), get_node(g,n2,i2), dir)

	end 

end

function GF_Decimation_fromGraph(Energy::Number, 
																 g::MetaDiGraph,
																 data_H::Dict,
																 translate::Function;
																 kwargs...
																 )::Function

	out1 = GF_Decimation_fromGraph(Energy, g, data_H; kwargs...)

	return function out2(args...; kwargs1...)::AbstractMatrix{ComplexF64}

		ni1ni2,slice = translate(unpack_argsG(args...)[1]...) 

		return view(out1(ni1ni2...; kwargs1...), slice...)

	end 

end 


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
	@assert !LA.ishermitian(H_intercell)  H_intercell

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

