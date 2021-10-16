module H_Superconductor

#############################################################################



import ..LA 

import ..Utils, ..TBmodel, ..Algebra




#===========================================================================#
#
# hopp
#
#---------------------------------------------------------------------------#

function ToBasis(Nambu::Bool=false,Zero=0)

  Nambu ? h->[h Zero; Zero -conj(h)] : h->h 

end

function BasisInfo(using_SC_basis::Bool)::Function

	!using_SC_basis && return (label=nothing;kwargs...) -> [1.0]

	function get_vector(label::AbstractString; kwargs...)::Vector{Float64}
#			kwargs:#		charge=false, spin=false, electron=false, hole=false)
	

			label=="Charge" && return get_vector(;charge=true)

			label=="QP" && return get_vector()

			label=="E" && return get_vector(;electron=true)

			label=="H" && return get_vector(;hole=true)

			label=="Sz" && return get_vector(;spin=true)

			error("Label $label not understood")

	end 


	function get_vector(label::Nothing=nothing; kwargs...)::Vector{Float64}

		out = map([(1,1), (1,-1), (-1,1), (-1,-1)]) do (C,S)

			mapreduce(*, (
										(:charge,		C),
										(:spin,			S),
										(:electron,	Int(C>0)),
										(:hole,			Int(C<0))
										)
								) do (k,v)

				get(kwargs, k, false) ? v : 1

			end 

		end 

		return length(unique(out))==1 ? out[1:1] : out 


	end 

end


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#



function SC_Domain(param_H_::NamedTuple, dist::AbstractVector{<:Real}; 
									 indomain=nothing, nr_uc::Int=1, 
									 dist_tol::Float64=1e-5, hopp_cutoff::Float64=1e-6,
									 )::Dict{Symbol,Any}


        # ---------------- initialize parameters --------------------- #      

  default_param_H = (
	Hopping			= 1.0,
	ChemicalPotential	= 2.0,
	LocalPotential 		= 0.0,
        SC_Gap			= nothing,
#	Electric_Field_x	= 0.0,
#	Electric_Field_y	= 0.0,
#	Electric_Field_z	= 0.0,
	Lattice_Imbalance	= 0.0,
#	Rashba_SOC		= 0.0,
					)


  param_H = Utils.Combine_NamedTuples(param_H_, default_param_H)


  using_SC_basis = !isnothing(param_H[:SC_Gap])

  PM = Algebra.PauliMatrices()

  samplehopp =  using_SC_basis ? PM[0] : rand(Complex{Float64},1,1)

  Zero,One = Base.zero(samplehopp),Base.one(samplehopp) 

  toBasis = ToBasis(using_SC_basis,Zero)

  d0 = size(toBasis(samplehopp),1)

  same = Utils.fSame(dist_tol)

        # ------------------- needed for the hoppings ----------------- #      


  Append,Sum,Hoppings = TBmodel.Add_Hopping_Terms(d0,hopp_cutoff)
		# value, condition, (function)
  
  maxdist,update = TBmodel.Update_Maxdist(dist_tol,hopp_cutoff)

  
        # -------------- nearest neighbor ------------------------ #      

  if length(dist)>0


    Hoppings = Append(Hoppings, param_H[:Hopping],
  		(ri,rj) -> same(ri.-rj,dist[1]),
  		toBasis(One),
  							)
    maxdist = update(maxdist,param_H[:Hopping],dist[1])

  end  

        # -------------- atom bias (lattice imbalance) ----------- #      

#  Hoppings = Append(Hoppings, param_H[:Lattice_Imbalance],
#		(ri,rj) -> same(ri,rj),
#		(ri,rj) -> 2*ri[sl] - sum(values(slcoord)),
#							)

        # ----------------- local potential --------------- #      

  if isa(param_H[:LocalPotential],Number)

    Hoppings = Append(Hoppings, param_H[:LocalPotential]
				+ param_H[:ChemicalPotential],
  		(ri,rj) -> same(ri,rj),
  		toBasis(One),
							)

  else

    Hoppings = Append(Hoppings, 1.0, 
  		(ri,rj) -> same(ri,rj),
  		(ri,rj) -> toBasis(One*(
				param_H[:LocalPotential](ri,rj)
				+ param_H[:ChemicalPotential]),
					),
							)
    
  end



        # -------------- superconducting pairing ----------------- #      


  if using_SC_basis

    val, nmax = param_H[:SC_Gap]

    Hoppings = Append(Hoppings, 1, true, val )

    nr_uc = max(Int64(nmax),nr_uc)

    maxdist = update(maxdist,1,nmax)


  end
        # ------------------- electric field -------------------- #      

#  Hoppings = Append(Hoppings, param_H[:Electric_Field_x],
#		(ri,rj) -> same(ri,rj),
#		(ri,rj) -> ri[1],
#							)
#
#
#  Hoppings = Append(Hoppings, param_H[:Electric_Field_y],
#		(ri,rj) -> same(ri,rj),
#		(ri,rj) -> ri[2],
#							)
#
#
#  Hoppings = Append(Hoppings, param_H[:Electric_Field_z]*(dim==3),
#		(ri,rj) -> same(ri,rj),
#		(ri,rj) -> ri[3],
#							)
#
#        # -------------- anti-Haldane nnn imaginary hopping ------ #      
#
#
#  nn(ri) = findall(same.(Algebra.FlatOuterDist(hcat(ri[1:2]...),AllAtoms[:,1:2]),dist[1]))
#
#  chirality(ri,rj) = sign(Algebra.cross(eachcol(hcat(ri,rj) .- AllAtoms[intersect(nn(ri),nn(rj))[1],:] )...)[3])
#
#
#  Hoppings = Append(Hoppings, param_H[:AntiHaldane],
#		(ri,rj) -> same(ri[1:2].-rj[1:2],dist[2]),
#		(ri,rj) -> 1im/sqrt(27)*ri[sl]*chirality(ri,rj),
#							)
#
#  maxdist = update(maxdist,param_H[:AntiHaldane],dist[2])

        # -------------- Kane-Mele spin-orbit coupling ----------- #      

#    if same(la.norm(dif[:2]),dist_nnn)*same(abs(dif[2]),0.)*neighbors.size == 0
#  push!(hoppings, TBmodel.Hopping_Term( 
#		(ri,rj) -> isapprox(LA.norm(ri.-rj),dist_nnn,atol=dist_tol),
#		param_H[:KaneMele_SOC],
#		,
#		tol = hopp_cutoff))

        # -------------- Rashba spin-orbit coupling -------------- #      
#
#  push!(hoppings, TBmodel.Hopping_Term( 
#		(ri,rj) -> isapprox(LA.norm(ri.-rj),dist_nnn,atol=dist_tol),
#		param_H[:Rashba_SOC],
#		,
#		tol = hopp_cutoff))
#


        # -------------- compute max_dist and sum up ------------- #      

  cond(ri,rj) = (isnothing(indomain) || indomain(ri,rj)) && LA.norm(ri.-rj) < maxdist
  
#  returnB_Input, (Sum(Hoppings,cond), hopp_cutoff)

	return Dict(
		:Hopping => Sum(Hoppings,cond),
		:Nr_UCs => nr_uc,	# provided the kind of hoppings above 
		:hopp_cutoff => hopp_cutoff,
		:dist_tol => dist_tol,
		:SC_basis => using_SC_basis,
		:BasisInfo => BasisInfo(using_SC_basis),
		:nr_orb => d0,
		)
			

end




#function ParticleHole_Flavours(Nambu)
#
#  samplehopp = rand(Complex{Float64},1+Nambu,1+Nambu) # add spin if Nambu
#
#  Zero,One = Base.zero(samplehopp),Base.one(samplehopp) 
#
#  return LA.diag(ToBasis(Nambu,Zero)(One)) # Nambu doubling
#
#end


#
#
#function ParticleHole_Operator(param_H,nr_at::Int64)
#
#  return TBmodel.Constant_Operator(nr_at,Charge_Operator_UC(param_H))
#end
#
#
#function Spin_Operator(param_H,atoms,axis=3)
#
#  if isnothing(param_H[:SC_Gap]) 
#
#    return TBmodel.Constant_Operator(size(atoms,1),hcat(0))
#
#  end
#
#  S = Algebra.PauliMatrices()[axis]
#
#  M = ToBasis(true,Base.zero(S))(S)
#
#  return TBmodel.Constant_Operator(size(atoms,1),M)
#
#end








###############################################################################

end

