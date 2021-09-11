module ObservablesFromGF



import ..LA

import ..GreensFunctions, ..Utils, ..Algebra, ..Operators, ..LayeredSystem


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
						 Op=[1], kwargs...)::Float64

	trace = Operators.Trace("orbitals", Op; sum_up=true, kwargs...)
	
	return trace(-1/pi*imag(LA.diag(Gr)))


end

#===========================================================================#
#
# Local DOS from a Green's function given as a function of layers
#
#---------------------------------------------------------------------------#


function LDOS_Decimation(GD::Function, NrLayers::Int, indsLayer::Function; 
												 Op=1, VirtLeads...)::Vector{Float64}

	dev_atoms = Dict(("Layer",L) => indsLayer(L) for L in 1:NrLayers)

	nr_at = sum(length, values(dev_atoms))

	lead_atoms = LayeredSystem.LeadAtomOrder(nr_at; VirtLeads...)

	ldos = zeros(Float64, mapreduce(length, +, values(lead_atoms), init=nr_at))


	for d in (dev_atoms,lead_atoms), (key,inds) in pairs(d)

		g = GD(key...)

		ldos[inds] = LDOS(g; Op=Op, nr_at=length(inds), size_H=size(g,1))

	end

	return ldos
	
end





function DOS_Decimation(GD::Function, NrLayers::Int, indsLayer::Function; 
												Op=[1], dim::Int, VirtLeads...)::Float64

	out = 0.0

	for d in (pairs(LayeredSystem.LeadAtomOrder(;VirtLeads...)),
						(("Layer",L)=>indsLayer(L) for L=1:NrLayers))

		for (key,inds) in d

			g = GD(key...) 

			out += DOS(GD(key...), Op=Op, nr_at=length(inds), size_H=size(g,1))
			
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

function JosephsonCurrent(GD::Function, i::Int; f=LA.tr)

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

function TunnelingConductance_LeeFisher(GD::Function, i::Int, j::Int=i;
																				f::Function=LA.tr)

	G(n,m) = (g -> (g' - g)/(2im))(GD(n,m))

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


function CaroliConductance(G1, source, drain, G2=G1; f=LA.tr)

#  sp(A) = SpA.dropzeros(SpA.sparse(A))
#
	GammaS = GreensFunctions.DecayWidth(source)
	
	GammaD = GreensFunctions.DecayWidth(drain)


#  out = f(LA.diag(sp(GammaS)*sp(G1)*sp(GammaD)*sp(G2)'))
##  out = LA.tr(GammaS*G1*GammaD*G2')
	
	return f(GammaS*G1*GammaD*G2')

#  Im = sum(abs.(imag(out)) ./ (abs.(out) .+1e-12))/length(out)

#  Im > 1e-6 && println("Caroli T has imaginary part. mean(imag/abs)=",round(Im,digits=7))

#  any(abs.(imag(out)) .> 1e-9) && println("T has imaginary part")


#  return real(out)

end




#	A = G*W*G'

#	GammaS = GreensFunctions.DecayWidth(source)
	
#	GammaD = GreensFunctions.DecayWidth(drain)

#	return f(GammaS*G1*GammaD*G2')


#h = H ( (layer2,atom2)-> (layer1,atom1))



function BondTij(Gi::AbstractMatrix{ComplexF64},
								 W::AbstractMatrix{ComplexF64},
								 Gj::AbstractMatrix{ComplexF64},
								 Hji::AbstractMatrix{ComplexF64},
								 f::Function)

	-2imag(f(	Gi*W*Gj'*Hji - Gj*W*Gi'*Hji'	))

end 



function BondTransmission(G::Function, 
													SE_lead, Bonds, RBonds, Hoppings; f=LA.tr, kwargs...)

	W = GreensFunctions.DecayWidth(SE_lead)

	bondT = zeros(length(Bonds))

	for sector in Utils.IdentifySectors(first.(sort(Bonds)))
							# for each atom, basically

		Gi = G( Bonds[sector[1]][1] )

		for bond_index in sector 

			j = Bonds[bond_index][2]

			bondT[bond_index] = BondTij(Gi, W, G(j), Hoppings[bond_index], f)

		end

	end 

	return bondT


																					
end


function addSiteTiTj!(siteT, (i,j), (Ri,Rj), Tij)

	t = Tij .* (Rj-Ri)

	siteT[ i, : ] .+=  t 
	siteT[ j, : ] .+=  t

end 

function SiteTransmission(BondT, Bonds, RBonds)

	siteT = zeros(Float64, maximum(maximum.(Bonds)), length(RBonds[1][1]) )

	for BRT in zip(Bonds, RBonds, BondT)

		addSiteTiTj!(siteT, BRT...)

	end

	return siteT 

end



function SiteTransmission(G, SE_lead, Bonds, RBonds, Hoppings; f=LA.tr, kwargs...)

	W = GreensFunctions.DecayWidth(SE_lead)


	siteT = zeros(Float64, maximum(maximum.(Bonds)), length(RBonds[1][1]) )

	for sector in Utils.IdentifySectors(first.(sort(Bonds)))
							# for each atom i, basically

		Gi = G( Bonds[sector[1]][1] )

		for bond_index in sector 
							# for each atom j connected to atom i

			(i,j) = Bonds[bond_index]

			addSiteTiTj!(siteT, (i,j), RBonds[bond_index],
									BondTij(Gi, W, G(j), Hoppings[bond_index], f),
									)
									 
		end


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
