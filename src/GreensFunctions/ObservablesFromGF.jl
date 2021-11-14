module ObservablesFromGF



import ..LA

import ..GreensFunctions, ..Utils, ..Algebra, ..Operators, ..LayeredSystem, ..ArrayOps


#===========================================================================#
#
# Local DOS from a Green's function given as a matrix
#
#---------------------------------------------------------------------------#


function LDOS(Gr::AbstractMatrix{ComplexF64}; 
							Op=1, kwargs...)::Vector{Float64}

	trace = Operators.Trace("orbitals", Op; sum_up=false, kwargs...)
	
	return trace(-1/pi*imag(LA.diag(Gr)))

end


function DOS(Gr::AbstractMatrix{ComplexF64}; 
						 Op=1, kwargs...)::Float64

	D = -1/pi*imag(LA.diag(Gr)) 

	length(Op)==1 && return only(Op)*sum(D) 

	return Operators.Trace("orbitals", Op; sum_up=true, kwargs...)(D)

end

#===========================================================================#
#
# Local DOS from a Green's function given as a function of layers
#
#---------------------------------------------------------------------------#


function LDOS_Decimation(GD::Function, NrLayers::Int, indsLayer::Function; 
												 Op=1, dim::Int, VirtLeads...)::Vector{Float64}


	dev_atoms = Dict(("Layer",L) => indsLayer(L) for L in 1:NrLayers)

	nr_at = sum(length, values(dev_atoms))

	lead_atoms = LayeredSystem.LeadAtomOrder(nr_at; dim=dim, VirtLeads...)

	ldos = zeros(Float64, mapreduce(length, +, values(lead_atoms), init=nr_at))


	for d in (dev_atoms,lead_atoms), (key,inds) in pairs(d)

		g = GD(key...)

		ldos[inds] = LDOS(g; Op=Op, nr_at=length(inds), size_H=size(g,1), dim=dim)

	end

	return ldos
	
end





function DOS_Decimation(GD::Function, NrLayers::Int, indsLayer::Function; 
												Op=[1], dim::Int, VirtLeads...)::Float64

	out = 0.0

	for d in (pairs(LayeredSystem.LeadAtomOrder(;dim=dim, VirtLeads...)),
						(("Layer",L)=>indsLayer(L) for L=1:NrLayers))

		for (key,inds) in d

			g = GD(key...) 

			out += DOS(g; Op=Op, nr_at=length(inds), size_H=size(g,1), dim=dim)
			
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
										G, NrLayers::Int, IndsAtomsOfLayer, Op=[1], VirtLeads=[];
										dos::Bool=true, ldos::Bool=true, 
										doskey=nothing, ldoskey=nothing,
										kwargs...)

	if ldos

		LDOS_ = LDOS_Decimation(G, NrLayers, IndsAtomsOfLayer; 
													 Op=Op, VirtLeads..., kwargs...)
		
		if !dos 

			isnothing(ldoskey) && return (nothing, LDOS_)

			return Dict(ldoskey=>LDOS_)

		else 	

			isnothing(doskey) | isnothing(ldoskey) && return (sum(LDOS_),LDOS_)

			return Dict(ldoskey=>LDOS_, doskey=>sum(LDOS_))

		end	

	elseif dos

		DOS = DOS_Decimation(G, NrLayers, IndsAtomsOfLayer; 
												 Op=Op, VirtLeads..., kwargs...)

		isnothing(doskey) && return (DOS, nothing)

		return Dict(doskey=>DOS)

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


function CaroliConductance(G::Function, get_SE::Function,
													 source, drain, uc::Int;
														f::Function=LA.tr
													 )::ComplexF64

	gSS = G(source, uc, source, uc)
	
	gSD = source==drain ? gSS : G(source, uc, drain, uc) 

	SigmaS, SigmaD = get_SE(source, uc), get_SE(drain, uc)

	A = real(f(Algebra.Commutator(gSS,SigmaS)))

	@assert isapprox(0, A, atol=1e-10) A

	out1 = CaroliConductance(gSD, SigmaS, SigmaD; f=f)

#	out2 = CaroliConductance2(gSD, SigmaS, SigmaD; f=f)

#	@assert isapprox(out1, out2, atol=1e-14) 

	return out1 

end 





function CaroliConductance(G1::AbstractMatrix, 
													 sigma_source::AbstractMatrix, 
													 sigma_drain::AbstractMatrix, 
													 G2::AbstractMatrix=G1; 
													 f::Function=LA.tr)::ComplexF64

	# formula (7) of Meir+Wingreen, PRL 68, 2512 (1992)  
	# or formula 
	#
	# of Caroli et al, J.Phys.C: Solid State Phys. 4 916 (1971)

	f(*(GreensFunctions.DecayWidth(sigma_source),
			G1,
			GreensFunctions.DecayWidth(sigma_drain),
			G2'
			)
		)

end



function CaroliConductance2(G::Function, get_SE::Function,
														source, drain, uc::Int;
														f::Function=LA.tr
#														kwargs...
														)::ComplexF64


	gSS = G(source, uc, source, uc)
	
	gSD = source==drain ? gSS : G(source, uc, drain, uc) 

	SigmaS, SigmaD = get_SE(source, uc), get_SE(drain, uc)


	A = real(f(Algebra.Commutator(gSS,SigmaS)))

	@assert isapprox(0, A, atol=1e-10) A


#	out1 = CaroliConductance(gSD, SigmaS, SigmaD; f=f)#kwargs...)

	out2 = CaroliConductance2(gSD, SigmaS, SigmaD; f=f)#kwargs...)


#	@assert isapprox(out1, out2, atol=1e-14)

	return out2

end  

function CaroliConductance2(Gr::Function, Ga::Function,
														get_SE::Function,
														source, drain, uc::Int;
														f::Function=LA.tr
#														kwargs...
														)::ComplexF64

	gSD_ret = Gr(source, uc, drain, uc) 

	gDS_adv = Ga(drain, uc, source, uc) 



	@assert isapprox(gDS_adv,gSD_ret') "NO TRS"

	SigmaS, SigmaD = get_SE(source, uc), get_SE(drain, uc)

#	out1 = CaroliConductance(gSD, SigmaS, SigmaD; f=f)

	out2 = CaroliConductance2(gSD_ret, SigmaS, SigmaD, gDS_adv; f=f)


#	@assert isapprox(out1, out2, atol=1e-14)

	return out2

end  




function CaroliConductance2(gSD_ret::AbstractMatrix,
														SigmaS::AbstractMatrix, 
														SigmaD::AbstractMatrix,
														gDS_adv...;
														f::Function=LA.tr)::ComplexF64 

	f(longitudinal_conductance_matrix(gSD_ret, SigmaS, SigmaD, gDS_adv...))


end 

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

	1im*Algebra.Commutator(SigmaS,
												 gSD_ret*GreensFunctions.DecayWidth(SigmaD)*gDS_adv,
												 4=>adjoint)

end 


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function addSiteTiTj!(siteT::T,
											dim::Int,
											inds::Tuple{Vararg{Int}}, 
											(Ri,Rj),
											Tij::Float64)::T where T<:AbstractMatrix{<:Real}

	t = Tij .* (Rj-Ri)


	for i in inds 

		selectdim(siteT, dim, i) .+= t 

	end 

	return siteT 

end 


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function BondTij(Gi::AbstractMatrix{ComplexF64},
								 W::AbstractMatrix{ComplexF64},
								 Gj::AbstractMatrix{ComplexF64},
								 Hji::AbstractMatrix{ComplexF64};
								 kwargs...
								)::Float64

#	-2.0*imag(f(Gi*W*Gj'*Hji - Gj*W*Gi'*Hji')) # wrong
# 
# "im" comes already from the difference (i->j) - (j->i)! 

	BondTij(Gi*W*Gj', Hji; kwargs...)

end 


function BondTij(GiW::AbstractMatrix{ComplexF64},
								 Gj::AbstractMatrix{ComplexF64},
								 Hji::AbstractMatrix{ComplexF64};
								 kwargs...)::Float64

	BondTij(GiW*Gj',Hji; kwargs...)

end 

function BondTij(Gi_W_Gjdagg::AbstractMatrix{ComplexF64},
								 Hji::AbstractMatrix{ComplexF64};
								 f::Function=LA.tr, kwargs...)::Float64

	-2.0*imag(f(Gi_W_Gjdagg*Hji))

end 




function BondTransmission(G::Function, 
													Hoppings::AbstractVector{<:AbstractMatrix},
													 Bonds::AbstractVector{NTuple{2,Int}},
													 SE_lead::Function, lead_args...;
													 kwargs...)::Vector{Float64}

	W = GreensFunctions.DecayWidth(SE_lead(lead_args...))

	bondT = zeros(Float64, length(Bonds))

	for sector in Utils.IdentifySectors(first.(Bonds))
							# for each atom, basically, provided Bonds are sorted

		i = Bonds[sector[1]][1]

		GiW = G("Atom", i, lead_args...)*W

		for bond_index in sector 

			j = Bonds[bond_index][2]

			bondT[bond_index] = BondTij(GiW,
																	G("Atom", j, lead_args...),
																	Hoppings[bond_index]; kwargs...)
		end

	end 

	return bondT

end




function SiteTransmission(G::Function, 
													Hoppings::AbstractVector{<:AbstractMatrix},
													 Bonds::AbstractVector{NTuple{2,Int}},
													 RBonds::AbstractVector{<:AbstractVector{<:AbstractVector{<:Real}}},
													 SE_lead::Function, lead_args...;
													dim::Int, kwargs...
													)::Matrix{Float64}

	siteT = ArrayOps.init_zeros(dim => maximum(maximum, Bonds),
															[2,1][dim] => length(RBonds[1][1]))


	W = GreensFunctions.DecayWidth(SE_lead(lead_args...))

	for sector in Utils.IdentifySectors(first.(Bonds))
							# for each atom, basically, provided Bonds are sorted

		i = Bonds[sector[1]][1]

		GiW = G("Atom", i, lead_args...)*W

		for bond_index in sector 

			j = Bonds[bond_index][2]

			bt = BondTij(GiW, G("Atom", j, lead_args...), Hoppings[bond_index];
									 kwargs...) 

			addSiteTiTj!(siteT, dim, (i,j), RBonds[bond_index], bt) 

		end

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
									f::Function=LA.tr)::Float64

 	-2real(f((G("Atom",i,"Atom",j) - G("Atom",j,"Atom",i)')*Hji))

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
													 dim::Int,
													 kwargs...)::Matrix{Float64}

	siteT = ArrayOps.init_zeros(dim => maximum(maximum, Bonds),
															[2,1][dim] => length(RBonds[1][1]))

	for (b,Rb,Hb) in zip(Bonds, RBonds, Hoppings)

		addSiteTiTj!(siteT, dim, b, Rb, BondTij0(G, b, Hb; kwargs...))

	end 

	return siteT 

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

	for BRT in zip(Bonds, RBonds, BondT)

		addSiteTiTj!(siteT, dim, BRT...)

	end

	return siteT 

end






#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#








































#############################################################################
end
