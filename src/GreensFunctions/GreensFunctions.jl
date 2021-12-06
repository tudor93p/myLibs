module GreensFunctions
#############################################################################

import ..LA, ..ArrayOps


import ..Utils, ..TBmodel, ..Lattices

import ..Graph, ..LayeredSystem


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

GF(g::AbstractMatrix) = g


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


function SelfEn_fromGDecim(G::Function, 
													 VirtLeads::AbstractDict, 
													 LeadLayerSlicer::Function
													 )::Function

	function self_en(k,i::Int=1)::Matrix

		(K,i),slice = LeadLayerSlicer(k,i)

		inter = VirtLeads[typeof(first(keys(VirtLeads)))(K)][:intercell]

		dir = if K=="LeftLead" 
			
			"left"
			
					elseif K=="RightLead"
						
			"right"

					end 


		return SelfEn(inter[min(i,end)]'[slice...,slice...],
																	G(string(k),i+1,dir=dir)
																		)
	end


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
		

		
	
	
	
	


	islead(n) = occursin("Lead",get_prop(:type,n))

	islayer(n) = get_prop(:type,n)=="Layer"


	
	
	function H(args...)

		ns = f(args...)
		
		(length(ns) == 1 || ns[1]==ns[2]) && return get_prop(:H,ns[1])
	
		Graph.has_prop(g,ns...,:H) && return get_prop(:H,ns...)
	
		rns = reverse(ns)
	
		Graph.has_prop(g,rns...,:H) && return get_prop(:H,rns...)'

		return zeros(size(get_prop(:H,ns[1]),1),size(get_prop(:H,ns[2]),1))

	end
	
	
	function setEnergy_LeadGF(g_,Energy)
	
		for ll in get_prop(:LeadLabels)
			
			gfs = get_prop(:GFf,f(ll,1)[1])(Energy)	
		
			for uc in 1:get_prop(:UCsLeads)

				Graph.set_prop!(g_,f(ll,uc)[1],:GF,gfs[min(uc,end)])
			end
		end
	
		return n -> Graph.get_prop(g_,n,:GF)
	end
	
	
	
	function left(n)
	
	  ns = Graph.in_neighbors(g,n)
	
		return isempty(ns) ? nothing : ns[1]
	
	end
	
	function right(n)
	
	  ns = Graph.out_neighbors(g,n)
	
		return isempty(ns) ? nothing : ns[1]
	
	end
	
	
	function next(n::Int64,dir::String) 
	
		dir == "right" && return right(n)
	
		dir == "left" && return left(n)
	
		error("dir should be left or right")
	
	
	end
	


	function meets_layer(n,dir::String="both")::Bool
	
		isnothing(n) && return false

		(dir == "both" || islayer(n)) && return true
	
		return meets_layer(next(n,dir),dir)
	
	end
	
	
	
	
	
	function bring_closer(n,m) 
	
		isempty(Graph.FindPath(g,n,m)) &&	return false,left(n),right(m)
	
		return true,right(n),left(m)
	
	end
	
	function lead_extends(n,dir)
	
		return islead(n) && (dir == "both" || !meets_layer(n,dir))
	
	end



	return 	setEnergy_LeadGF,
					islead,
					meets_layer,
					lead_extends,
					H,
					next,
					bring_closer,
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

	HoppMatr(args...) = TBmodel.HoppingMatrix(AtomsOfLayer.(args)...;Hopping...)

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

	@assert haskey(kwargs,:leads_have_imag)
	leads_have_imag = kwargs[:leads_have_imag]





	setEnergy_LeadGF,
	islead,
	meets_layer,
	lead_extends,
	H,
	next,
	bring_closer,
	node	= GraphLayeredSystem_Utils(g)


	reEnergy = leads_have_imag ? Energy::Real : real(Energy::Complex)


	SemiInfLeadGF = setEnergy_LeadGF(g,Energy)


	# ----- the dictionary where the values are stored -------- #
	
	Gd = Dict{Tuple{Int64,Int64,String},Array{Complex{Float64},2}}()



	function G(n::Int,m::Int,dir::AbstractString)::AbstractMatrix
	
		n==m && islead(n) && !meets_layer(n,dir) && return SemiInfLeadGF(n)

							#	this should not be stored, already stored in the graph!

		n==m && return get!(Gd,(n,m,dir)) do 
																			Gnn(n,dir) end

		return get!(Gd,(n,m,dir)) do 
																Gnm(n,m,dir) end

	end


	# ----- methods to calculate the GF ------- #	

	function Gnn(n::Int, dir::AbstractString)::Matrix{ComplexF64}
	
		lead_extends(n,dir) && return GF(SemiInfLeadGF(n),
																			coupling_toSystem(n,dir)...)
	
		return GF(reEnergy, H(n), coupling(n,dir)... ) 
	
	end

	function Gnm(n::Int, m::Int, dir::AbstractString)::Matrix{ComplexF64}

		isnleft, n1, m1  = bring_closer(n,m)
	
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
	
		for dir in ["right","left"]
	
			meets_layer(src,dir) && return coupling(src,dir)
	
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
									kwargs...)::AbstractMatrix

			out1(name1, index1, name2, index2)

		end 

		function out1(ni1ni2::Tuple{Tuple,Tuple}; kwargs...)::AbstractMatrix

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
									kwargs...)::AbstractMatrix

			out2(name1, index1, name2, index2)

		end 

		function out2(ni1ni2::Tuple{Tuple,Tuple}; kwargs...)::AbstractMatrix

			out2(ni1ni2...)

		end

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

	first ∘ GF_Surface(latt, [target_arg], args...)


end

function GF_Surface(latt::Lattices.Lattice,
										target_arg::AbstractVector{<:AbstractString},
										hopping::AbstractDict,
										delta::Real,
										)::Function

#	gf_ = GF_Surface(latt, target_arg, hopping)

	intra, inter = TBmodel.BlochHamilt_ParallPerp(latt, Lattices.NearbyUCs, hopping)
	
	@assert Lattices.LattDim(latt)==1 

	intra_, inter_ = intra(), only(values(inter))
		

	return function gf(E::Real)::Vector{Matrix{ComplexF64}}
		
		gfs = GreensFunctions.GF_SanchoRubio(E+delta*im, intra_, inter_;
																				 target=join(target_arg))

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

	@assert LA.ishermitian(H_intracell)  

	Errors = zeros(max_iter)


	alpha = copy(Matrix(H_intercell))
	beta = copy(alpha')


	eps = cat((H_intracell for i in 0:(e2|e3))..., dims=3)

	gBulk = GF(Energy, selectdim(eps, 3, 1))


	for iter in 1:max_iter  

		SanchoRubio_converged!(Errors, iter, alpha, beta; kwargs...) && break  

		ga = gBulk*alpha  

		gb = gBulk*beta
		
		selectdim(eps, 3, 1) .+= beta*ga 

		eps .+= (alpha * gb)[:,:,:] 

		alpha *= ga  

		beta *= gb 

		gBulk = GF(Energy, selectdim(eps, 3, 1))

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

function PrepareLead(atoms::AbstractMatrix{<:Number})::Dict{Symbol,Any}

	Dict{Symbol,Any}(:head=>[atoms])

end 


function PrepareLead(label::AbstractString, 
										 LeadAtoms::AbstractMatrix{<:Number}
										 )::Dict{Symbol,Any}

	merge!(Dict{Symbol,Any}(:label => label), PrepareLead(LeadAtoms))
	
end  


function PrepareLead(label::AbstractString, Lead::Lattices.Lattice
										 )::Dict{Symbol,Any}

	PrepareLead(label, Lattices.PosAtoms(Lead))

end 


function PrepareLead(label::AbstractString, 
										 Lead::Union{Lattices.Lattice, AbstractMatrix{<:Number}},
										 BridgeAtoms::AbstractMatrix{<:Number}
										 )::Dict{Symbol,Any}

	merge!(vcat, PrepareLead(BridgeAtoms), PrepareLead(label, Lead)) 

end 

function PrepareLead(label::AbstractString,
										 lead_latt::Lattices.Lattice,
										 lead_hopp::AbstractDict,
										 delta::Vararg{Float64},
										 )::Dict{Symbol,Any} 

	PrepareLead(label, lead_latt, lead_hopp, lead_hopp, delta...)

end 

function PrepareLead(label::AbstractString,
										 lead_latt::Lattices.Lattice,
										 coupl_hopp::AbstractDict,
										 lead_hopp::AbstractDict,
										 delta::Vararg{Float64},
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


	return merge!(PrepareLead(label, LeadAtoms0), Dict(	

						:coupling => coupling,

						:intracell => [LeadHoppMatr(LeadAtoms0)],

						:intercell => [LeadHoppMatr(LeadAtoms0, LeadAtoms1)],

						:GF => function GF(E::Number)::Vector{Matrix} 
						
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
				
				
	out[:GF] = function GF(E::Number)::Vector{Matrix}
	
				g = LeadGF(E)

				g1 = GF(E, BridgeIntra, BridgeToLead', g)
				
				g11 = inv(Array(E*LA.I - BridgeIntra - BridgeToLead*g*BridgeToLead'))

				@show g1≈g11
				@assert g1≈g11 LA.norm(g1-g11)

				return [g1, g]

		end

	return out

end



























































#############################################################################
end

