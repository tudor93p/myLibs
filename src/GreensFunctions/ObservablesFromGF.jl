module ObservablesFromGF



import ..LA

import ..GreensFunctions, ..Utils, ..Algebra, ..Operators, ..LayeredSystem, ..ArrayOps


#===========================================================================#
#
# Local DOS from a Green's function given as a matrix
#
#---------------------------------------------------------------------------#


function LDOS(Gr::AbstractMatrix{ComplexF64}; 
							proj::Function=identity, kwargs...)::Vector{Float64}

	trace = Operators.Trace(:orbitals; sum_up=false, kwargs...)
	
	return trace(proj(-1/pi*imag(LA.diag(Gr))))

end


function DOS(Gr::AbstractMatrix{ComplexF64}; 
						 proj::Function=identity, kwargs...)::Float64

	sum(proj(-1/pi*imag(LA.diag(Gr))))

#	length(Op)==1 && return only(Op)*sum(D) 

#	return Operators.Trace("orbitals", Op; sum_up=true, kwargs...)(D)

end

#===========================================================================#
#
# Local DOS from a Green's function given as a function of layers
#
#---------------------------------------------------------------------------#

function LDOS_Decimation(GD::Function; 
												 NrLayers::Int, 
												 IndsAtomsOfLayer::Function,
												 kwargs...)::Vector{Float64}

	LDOS_Decimation(GD, NrLayers, IndsAtomsOfLayer; kwargs...)

end

function LDOS_Decimation(GD::Function, NrLayers::Int, indsLayer::Function; 
												 proj::Function=identity, dim::Int, VirtLeads...)::Vector{Float64}


	dev_atoms = Dict(("Layer",L) => indsLayer(L) for L in 1:NrLayers)

	nr_at = sum(length, values(dev_atoms))

	lead_atoms = LayeredSystem.LeadAtomOrder(nr_at; dim=dim, VirtLeads...)

	ldos = zeros(Float64, mapreduce(length, +, values(lead_atoms), init=nr_at))


	for d in (dev_atoms,lead_atoms), (key,inds) in pairs(d)

		g = GD(key...)

		ldos[inds] = LDOS(g; proj=proj, nr_at=length(inds), size_H=size(g,1), dim=dim)

	end

	return ldos
	
end



function DOS_Decimation(GD::Function; 
												 NrLayers::Int, 
												 IndsAtomsOfLayer::Function,
												 kwargs...)::Float64

	DOS_Decimation(GD, NrLayers, IndsAtomsOfLayer; kwargs...)

end


function DOS_Decimation(GD::Function, NrLayers::Int, indsLayer::Function; 
												proj::Function=identity, dim::Int, VirtLeads...)::Float64

	out = 0.0

	for d in (pairs(LayeredSystem.LeadAtomOrder(;dim=dim, VirtLeads...)),
						(("Layer",L)=>indsLayer(L) for L=1:NrLayers))

		for (key,inds) in d

			g = GD(key...) 

			out += DOS(g; proj=proj, nr_at=length(inds), size_H=size(g,1), dim=dim)
			
		end 

	end


	return out

end




#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function ComputeDOSLDOS_Decimation(
										G::Function, NrLayers::Int, IndsAtomsOfLayer::Function, 
										proj::Function=identity, VirtLeads=[];
										dos::Bool=true, ldos::Bool=true, 
										doskey=nothing, ldoskey=nothing,
										kwargs...)

	if ldos

		LDOS_ = LDOS_Decimation(G, NrLayers, IndsAtomsOfLayer; 
													 proj=proj, VirtLeads..., kwargs...)
		
		if !dos 

			isnothing(ldoskey) && return (nothing, LDOS_)

			return Dict(ldoskey=>LDOS_)

		else 	

			isnothing(doskey) | isnothing(ldoskey) && return (sum(LDOS_),LDOS_)

			return Dict(ldoskey=>LDOS_, doskey=>sum(LDOS_))

		end	

	elseif dos

		DOS_ = DOS_Decimation(G, NrLayers, IndsAtomsOfLayer; 
												 proj=proj, VirtLeads..., kwargs...)

		isnothing(doskey) && return (DOS_, nothing)

		return Dict(doskey=>DOS_)

	end

	!isnothing(doskey) | !isnothing(ldoskey) && return Dict()

	return (nothing,nothing)

end



#===========================================================================#
#
# Josephson current from Furusaki, Physica B 203 (1994) 214-218 -- Eq. (5)
# 		see also Asano, PRB 63, 052512 (2001) -- Eq. (20)
#				PRB 74, 064507 (2006) -- Eq. (27) 
#
#---------------------------------------------------------------------------#

function JosephsonCurrent(GD::Function, i::Int; f::Function=LA.tr)

	@warn "f kwarg"

  -1im*f(GD(i,i-1) - GD(i-1,i))

  # GD = Green's function with Matsubara frequency wn
  # Furusaki:-- "Although the summation over wn in Eq. (5) i
  #		is originally from n=-oo to n=oo, we can alternatively sum up
  #		over positive wn only and take the real part of the sum. 
  #		Thus, we assume wn > 0 in the following discussion."
  # 	     -- "(...) calculate the current (...) in the normal region,
  #		[where] the electric charge is always conserved. 
  #		If the current is calculated in a superconducting region, 
  #		the source term proportional to the order parameter 
  #		must be included."

#  Im = sum(abs.(imag(out)) ./ (abs.(out) .+1e-12))/length(out)

#  Im > 1e-6 && println("J has imaginary part. mean(imag/abs)=",round(Im,digits=7))

# return vcat(real(out)...)
end

#===========================================================================#
#
# Tunneling conductance according to PRB 86, 174512 (2012),
#	formula based on Lee-Fisher: PRL 47, 882 (1981)
#
#---------------------------------------------------------------------------#

function TunnelingConductance_LeeFisher(GD::Function, 
																				lead, 
																				i::Int=2, j::Int=i;
																				f::Function=LA.tr)::ComplexF64

	@warn "f kwarg"

	function G(n::Int,m::Int)::AbstractMatrix{ComplexF64}
																					
		g = GD(lead, n, lead, m)

		return (g' - g)/(2im)

	end

  return f(	G(i,j)*G(j-1,i-1) 
									+ G(i-1,j-1)*G(j,i) 
									- G(i,j-1)*G(j,i-1) 
									- G(i-1,j)*G(j-1,i)
					)
	# prefactor aside, Lee-Fisher formula (3) with j0->i, j0'->j
	# becomes formula (15) in PRB with i=j=x+1


#	Im = sum(abs.(imag(La.diag(out))) ./ (abs.(LA.diag(out)) .+1e-12))/length(LA.diag(out))

#  Im > 1e-6 && println("LeeFisher T has imaginary part. mean(imag/abs)=",round(Im,digits=7))

#  return vcat(real(out)...)

end

#===========================================================================#
#
# Quantum conductance
#
#---------------------------------------------------------------------------#

function test_CaroliConductance(gSS::AbstractMatrix,
																SigmaS::AbstractMatrix;
																#f::Function=LA.tr
																)

	A = real(LA.tr(Algebra.Commutator(gSS, SigmaS)))

	@assert isapprox(0, A, atol=1e-10) A

end 




function CaroliConductance(G::Function, get_SE::Function,
													 source::AbstractString, drain::AbstractString, 
													 projectors::Vararg{Function};
													 kwargs...
													 )::Float64 

	CaroliConductance(G, get_SE, source, drain, 1, projectors...; kwargs...)

end 


function CaroliConductance(G::Function, get_SE::Function,
													 source::AbstractString, drain::AbstractString, 
													 uc::Int,
													 projectors::Vararg{Function};
													 kwargs...
													 )::Float64
	
	source==drain && return CaroliConductance(G, get_SE, source, uc, 
																						projectors...; kwargs...)

	
	SigmaS = get_SE(source, uc)

	test_CaroliConductance(G(source, uc, source, uc), SigmaS; kwargs...)


	return CaroliConductance(G(source, uc, drain, uc),
													 GreensFunctions.DecayWidth(SigmaS),
													 GreensFunctions.DecayWidth(get_SE(drain, uc)),
													 projectors...;
													 kwargs...)
end 



function CaroliConductance(G::Function, get_SE::Function,
													 lead::AbstractString, 
													 projectors::Vararg{Function};
													 kwargs...
													 )::Float64

	CaroliConductance(G, get_SE, lead, 1, projectors...; kwargs...)

end 


function CaroliConductance(G::Function, get_SE::Function,
													 lead::AbstractString, uc::Int,
													 projectors::Vararg{Function};
													 kwargs...
													 )::Float64

	g = G(lead, uc) 
	
	Sigma = get_SE(lead, uc)

	test_CaroliConductance(g, Sigma; kwargs...) 

	return CaroliConductance(g, 
													 GreensFunctions.DecayWidth(Sigma),
													 projectors...; kwargs...) 

end  

function CaroliConductance(G::AbstractMatrix, 
													 gamma::AbstractMatrix, 
													 projectors::Vararg{Function};
													 kwargs...)::Float64 

	CaroliConductance(G, gamma, gamma, projectors...; kwargs...)

end 


function CaroliConductance(Gr::Function, Ga::Function,
														get_SE::Function,
														source::AbstractString, drain::AbstractString, 
													 projectors::Vararg{Function};
														kwargs...
														)::ComplexF64

	CaroliConductance(Gr, Ga, get_SE, source, drain, 1,
										projectors...; kwargs...)
end 

function CaroliConductance(Gr::Function, Ga::Function,
														get_SE::Function,
														source::AbstractString, drain::AbstractString, 
														uc::Int,
													 projectors::Vararg{Function};
														kwargs...
														)::ComplexF64

	source==drain && return CaroliConductance(Gr, get_SE,
																						source, uc, projectors...;
																						kwargs...)

	gSD_ret = Gr(source, uc, drain, uc) 

	gDS_adv = Ga(drain, uc, source, uc) 

	SigmaS = get_SE(source, uc)


	#	@assert isapprox(gDS_adv,gSD_ret',rtol=1e-5) "NO TRS"
	
	isapprox(gDS_adv,gSD_ret',rtol=1e-5) || @warn "NO TRS?"

	test_CaroliConductance(Gr(source, uc, source, uc), SigmaS; kwargs...)



	return CaroliConductance(gSD_ret,
													 GreensFunctions.DecayWidth(SigmaS),
													 GreensFunctions.DecayWidth(get_SE(drain, uc)),
													 gDS_adv,
													 projectors...; kwargs...)

end  





function CaroliConductance(gR_SD::AbstractMatrix, 
													 gamma_source::AbstractMatrix, 
													 gamma_drain::AbstractMatrix, 
													 gA_DS::AbstractMatrix; 
#													 f::Function=LA.tr
													 )::ComplexF64

	@assert LA.ishermitian(gamma_source) && LA.ishermitian(gamma_drain)

	# formula (7) of Meir+Wingreen, PRL 68, 2512 (1992)  
	# or formula 
	#
	# of Caroli et al, J.Phys.C: Solid State Phys. 4 916 (1971)

	LA.tr(gamma_source * gR_SD * gamma_drain * gA_DS) 


end

function CaroliConductance(G1::AbstractMatrix, 
													 gamma_source::AbstractMatrix, 
													 gamma_drain::AbstractMatrix, 
													 G2::AbstractMatrix,
													 proj::Function, 
													 kwargs...)::ComplexF64  
	
	CaroliConductance(G1, gamma_source, gamma_drain, G2, proj, proj; kwargs...)

end

function CaroliConductance(G1::AbstractMatrix, 
													 gamma_source::AbstractMatrix, 
													 gamma_drain::AbstractMatrix, 
													 G2::AbstractMatrix,
													 proj_source::Function, proj_drain::Function;
													 kwargs...)::ComplexF64 

	CaroliConductance(G1, 
										proj_source(gamma_source),
										proj_drain(gamma_drain),
										G2;
										kwargs...)

end 



function CaroliConductance(G::AbstractMatrix, 
													 gamma_source::AbstractMatrix, 
													 gamma_drain::AbstractMatrix,
													 projectors::Vararg{Function};
													 kwargs...)::Float64

	real(CaroliConductance(G, gamma_source, gamma_drain, G',
												 projectors...; kwargs...))

end 




#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#










function FanoFactor2(G::Function, get_SE::Function,
														source, drain, uc::Int;
														kwargs...
														)::ComplexF64 

	FanoFactor2(G(source, uc, drain, uc), 
							get_SE(source, uc), 
							get_SE(drain, uc);
							kwargs...)

end

function FanoFactor2(args...; f::Function=LA.tr)::ComplexF64 

	X = longitudinal_conductance_matrix(args...)

	return 1 - f(X*X)/f(X)

end


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


#function LongitudinalBondTransmission(gSD::AbstractMatrix,
#																		 SigmaS::AbstractMatrix,
#																		 SigmaD::AbstractMatrix;
#																			)

#	[(i,i+1) for 
#	(1,2)
#	(2,3)
#	...
#	(nr_at-1,nr_at)

#	tr_orb = if haskey(kwargs,:f) 
#
#							kwargs[:f]::Function 
#
#					elseif haskey(kwargs, :nr_orb) 
#
#							kwargs[:nr_orb]::Int 
#
#							Operators.Trace(:orbitals; size_H=size(SigmaS,1), kwargs...)
#
#					else error() end 
#



#end  





function longitudinal_conductance_matrix(gSD_ret::AbstractMatrix,
																				 SigmaS::AbstractMatrix,
																				 SigmaD::AbstractMatrix,
																				 gDS_adv::AbstractMatrix=gSD_ret',
																				 )::Matrix{ComplexF64}

#	left: source, SigmaS, GammaS
#	right: drain, SigmaD, GammaD 
	# formula (10) of PRB 68, 075306 (2003) -- or (A21)

#	@assert LA.ishermitian(gamma_source) && LA.ishermitian(gamma_drain)


	1im*Algebra.Commutator(SigmaS,
												 gSD_ret*GreensFunctions.DecayWidth(SigmaD)*gDS_adv,
												 4=>adjoint)

end 


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function addSiteTiTj!(siteT::AbstractMatrix{<:Real},
											dim::Int,
											inds::NTuple{N,Int} where N,
											(Ri,Rj),
											Tij::Float64)::Nothing 

	t = Tij .* (Rj-Ri)

	for i in inds 

		selectdim(siteT, dim, i) .+= t 

	end 

	return 

end 


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function BondTij(Gr_iL::AbstractMatrix{ComplexF64},
								 W_L::AbstractMatrix{ComplexF64},
								 Ga_Lj::AbstractMatrix{ComplexF64},
								 Hji::AbstractMatrix{ComplexF64};
								 kwargs...
								)::Float64

#	-2.0*imag(f(Gi*W*Gj'*Hji - Gj*W*Gi'*Hji')) # wrong
# 
# "im" comes already from the difference (i->j) - (j->i)! 

	BondTij(Gr_iL * W_L * Ga_Lj, Hji; kwargs...)

end 


function BondTij(GriW::AbstractMatrix{ComplexF64},
								 Gaj::AbstractMatrix{ComplexF64},
								 Hji::AbstractMatrix{ComplexF64};
								 kwargs...)::Float64

	BondTij(GriW*Gaj,Hji; kwargs...)

end 

function BondTij(Gri_W_Gaj::AbstractMatrix{ComplexF64},
								 Hji::AbstractMatrix{ComplexF64};
#								 f::Function=LA.tr, 
								 kwargs...)::Float64

	-2.0*imag(LA.tr(Gri_W_Gaj * Hji))

end 





function BondTransmission_(
													Hoppings::AbstractVector{<:AbstractMatrix},
													 Bonds::AbstractVector{NTuple{2,Int}},
													 SE_lead::Function, 
													 lead::AbstractString,
													 lead_uc::Int,
													 Gs::Vararg{Function}
													 ;
													 kwargs...)::Vector{Float64}

	W = GreensFunctions.DecayWidth(SE_lead(lead,lead_uc))

	bondT = zeros(Float64, length(Bonds))

	for sector in Utils.IdentifySectors(first, Bonds)
							# for each atom, basically, provided Bonds are sorted

		BondTiJ!(sector, bondT, Hoppings, Bonds, Gs..., W, 
						 lead, lead_uc; kwargs...)

	end 

	return bondT

end


function BondTiJ!(sector::AbstractVector{Int},
									bondT::AbstractVector{Float64},
									Hoppings::AbstractVector{<:AbstractMatrix},
									Bonds::AbstractVector{NTuple{2,Int}},
									args...; kwargs...)::Nothing 

	BondTiJ!(view(bondT,sector), view(Hoppings, sector), view(Bonds,sector),
					 Bonds[sector[1]][1], args...; kwargs...)

end  

function BondTiJ(sector::AbstractVector{Int},
								#	bondT::AbstractVector{Float64},
									Hoppings::AbstractVector{<:AbstractMatrix},
									Bonds::AbstractVector{NTuple{2,Int}},
									args...; kwargs...)::Vector{Float64}

	bondT = zeros(length(sector))

	BondTiJ!(bondT, view(Hoppings, sector), view(Bonds,sector),
					 Bonds[sector[1]][1], args...; kwargs...)

	return bondT

end  


function BondTiJ!(bondT::AbstractVector{Float64},
									Hoppings::AbstractVector{<:AbstractMatrix},
									Bonds::AbstractVector{NTuple{2,Int}},
									i::Int,
									Gr::Function, #Ga::Function,
									SE_lead::AbstractMatrix,
									lead...; kwargs...)::Nothing
	
	BondTiJ!(bondT, Hoppings, Bonds, Gr("Atom",i,lead...)*SE_lead, 
					 adjoint∘Gr∘reverse, lead...; kwargs...)

end 

function BondTiJ!(bondT::AbstractVector{Float64},
									Hoppings::AbstractVector{<:AbstractMatrix},
									Bonds::AbstractVector{NTuple{2,Int}},
									i::Int,
									Gr::Function, Ga::Function,
									SE_lead::AbstractMatrix,
									lead...; kwargs...)::Nothing
	
	gr = Gr("Atom",i,lead...)
	
	@assert gr ≈ Ga(lead..., "Atom", i)'

	BondTiJ!(bondT, Hoppings, Bonds, gr*SE_lead, Ga, lead...; kwargs...)

end 



function BondTiJ!(bondT::AbstractVector{Float64},
									Hoppings::AbstractVector{<:AbstractMatrix},
									Bonds::AbstractVector{NTuple{2,Int}},
									GiW::AbstractMatrix{ComplexF64},
									Ga::Function,
									lead...; kwargs...)::Nothing

	for (ibt, (H,(i,j))) in enumerate(zip(Hoppings,Bonds))

		bondT[ibt] = BondTij(GiW, Ga((lead,("Atom",j))), H; kwargs...)

	end 

	return 

end  



function BondTransmission(Gr::Function, 
													Hoppings::AbstractVector{<:AbstractMatrix},
													args...; kwargs...)::Vector{Float64}

	BondTransmission_(Hoppings, args..., Gr; kwargs...)

end



function BondTransmission(Gr::Function, Ga::Function,
													args...; kwargs...)::Vector{Float64}

	BondTransmission_(args..., Gr, Ga; kwargs...)

end



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#




function SiteTransmission(G::Function, 
													Hoppings::AbstractVector{<:AbstractMatrix},
													args...; kwargs...)::Matrix{Float64}

	SiteTransmission_(Hoppings, args..., G; kwargs...)

end
 



function addSiteTiTJ!(sector::AbstractVector{Int},
											siteT::AbstractMatrix{Float64},
											Hoppings::AbstractVector{<:AbstractMatrix},
											Bonds::AbstractVector{NTuple{2,Int}},
											RBonds::AbstractVector{<:AbstractVector{<:AbstractVector{<:Real}}},
											args...; kwargs...)::Nothing

	SiteTransmission_fromBondT!(siteT,
															view(Bonds,sector), view(RBonds, sector), 
															BondTiJ(sector, Hoppings, Bonds, args...; 
																			kwargs...); kwargs...)

end 

function SiteTransmission(Gr::Function, Ga::Function,
													Hoppings::AbstractVector{<:AbstractMatrix},
													args...; kwargs...
													)::Matrix{Float64}

	SiteTransmission_(Hoppings, args..., Gr, Ga; kwargs...)

end 

function SiteTransmission_(Hoppings::AbstractVector{<:AbstractMatrix},
													 Bonds::AbstractVector{NTuple{2,Int}},
													 RBonds::AbstractVector{<:AbstractVector{<:AbstractVector{<:Real}}},
													 SE_lead::Function, 
													 lead::AbstractString,
													 lead_uc::Int,
													 Gs::Vararg{Function}
													 ;
													 dim::Int,
													 kwargs...)::Matrix{Float64}

	siteT = ArrayOps.init_zeros(dim => maximum(maximum, Bonds),
															[2,1][dim] => length(RBonds[1][1]))

	W = GreensFunctions.DecayWidth(SE_lead(lead, lead_uc))

	for sector in Utils.IdentifySectors(first, Bonds)

		addSiteTiTJ!(sector, siteT, Hoppings, Bonds, RBonds, Gs..., W, 
								 lead, lead_uc; dim=dim, kwargs...) 

	end 

	return siteT

end




#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function BondTij0(G::Function,
									(i,j)::NTuple{2,Int},
									Hji::AbstractMatrix{ComplexF64};
#									f::Function=LA.tr,
									kwargs...)::Float64

 	-2real(LA.tr((G("Atom",i,"Atom",j) - G("Atom",j,"Atom",i)')*Hji))

end 


function BondTransmission0(G::Function, 
													 Hoppings::AbstractVector{<:AbstractMatrix},
													 Bonds::AbstractVector{NTuple{2,Int}};
													 kwargs...)::Vector{Float64}

	bondT = zeros(Float64, length(Bonds))
	
	for (bond_index, bH) in enumerate(zip(Bonds,Hoppings))

		bondT[bond_index] = BondTij0(G, bH...; kwargs...) 

	end 

	return bondT 

end



function SiteTransmission0(G::Function, 
													 Hoppings::AbstractVector{<:AbstractMatrix},
													 Bonds::AbstractVector{NTuple{2,Int}},
													 RBonds::AbstractVector{<:AbstractVector{<:AbstractVector{<:Real}}};
													 kwargs...)::Matrix{Float64}

	SiteTransmission_fromBondT(BondTransmission0(G, Hoppings, Bonds; kwargs...),
														 Bonds, RBonds; kwargs...)


end














#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#






function SiteTransmission_fromBondT(BondT::AbstractVector{<:Real},
													 Bonds::AbstractVector{NTuple{2,Int}},
													 RBonds::AbstractVector{<:AbstractVector{<:AbstractVector{<:Real}}};
													 dim::Int)::Matrix{Float64}

	siteT = ArrayOps.init_zeros(dim => maximum(maximum, Bonds),
															[2,1][dim] => length(RBonds[1][1]))

	SiteTransmission_fromBondT!(siteT, Bonds, RBonds, BondT; dim=dim)
															
	return siteT 

end

function SiteTransmission_fromBondT!(siteT::AbstractMatrix{Float64},
																		 args...; dim::Int)::Nothing 

	for BRT in zip(args...)

		addSiteTiTj!(siteT, dim, BRT...)

	end

	return 

end






#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#








































#############################################################################
end
