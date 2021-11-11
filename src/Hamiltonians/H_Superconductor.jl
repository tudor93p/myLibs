module H_Superconductor

#############################################################################



import ..LA 

import ..Utils, ..TBmodel, ..Algebra




#===========================================================================#
#
# hopp
#
#---------------------------------------------------------------------------#

function ToBasis(Nambu::Bool=false, Z...)::Function 

	#  Nambu ? h->[h Zero; Zero -conj(h)] : h->h  
	
	ToBasis(Val(Nambu), Z...)


end

function ToBasis(::Val{false}, args...)::Function 

	TBmodel.matrixval 

end  

function ToBasis(::Val{true}, Zero::Tz=0.0im
								)::Function where Tz<:Union{Number, AbstractMatrix}

	function tobasis(h::Union{Number,AbstractMatrix})::Matrix{ComplexF64}

		[h Zero; Zero -conj(h)]
		
	end 
	
	tobasis(h::Function)::Function = tobasis ∘ h


	function tobasis(v::Union{Number, AbstractMatrix}, 
									 h::Union{Number,AbstractMatrix})::Matrix{ComplexF64}

		tobasis(TBmodel.matrixval(v,h))

	end 

	function tobasis(v::Union{Number, AbstractMatrix}, 
									 h::Function)::Function 

		tobasis∘TBmodel.matrixval(v, h)

	end 


	return tobasis 

	
end 





#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#



function BasisInfo(using_SC_basis::Bool)::Function

	BasisInfo(Val(using_SC_basis))
	
end 


function BasisInfo(using_SC_basis::Val{false})::Function

	function One(label=nothing; kwargs...)::Vector{Float64} 
		
		[1.0]

	end 

end 

function BasisInfo(using_SC_basis::Val{true})::Function

	kw0 = Dict(
						"Charge" => :charge,

						"QP" =>  nothing,

						"E" => :electron,

						"H" => :hole, 

						"Sz" => :spin

						)


	function get_vector(label::AbstractString; kwargs...)::Vector{Float64}

		@assert haskey(kw0,label) "Label '$label' not understood" 

		K = kw0[label]

		return get_vector(; (isnothing(K) ? () : Dict(K=>true))...)

	end 


	function get_vector(; kwargs...)::Vector{Float64}

		out = map([(1,1), (1,-1), (-1,1), (-1,-1)]) do (C,S)

			mapreduce(*, (
										(:charge,		C),
										(:spin,			S),
										(:electron,	Int(C>0)),
										(:hole,			Int(C<0))
										)
								) do (k,v)

				use = get(kwargs, k, false) 
			
				isa(use,Bool) && use && return v 
				
				isa(use,Int) && use==1 && return v 

				return 1 

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

function args_local_potential(LocalPotential::Number, 
						 ChemicalPotential::Number, 
						 )::Tuple{Float64, Float64}

	(LocalPotential + ChemicalPotential, 1.0)

end 



function args_local_potential(LocalPotential::Function, 
															ChemicalPotential::Number,
															)::Tuple{Float64, Function}

	(1.0, 

	 function total_local_pot(ri::AbstractVector, rj::AbstractVector
																	 )::AbstractMatrix
		 
		 LocalPotential(ri,rj) + ChemicalPotential
	 
	 end 

	 )

end 

function args_local_potential(;LocalPotential::Union{Number,Function},
															 ChemicalPotential::Number,
															 kwargs...)

	args_local_potential(LocalPotential, ChemicalPotential)

end 


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function PeierlsIntegral(A::AbstractVector, B::AbstractVector,
												 ::Val{:x})::Float64

	# Landau gauge A = [-y, 0, 0] 
	

		(A[1]-B[1])*(A[2]+B[2])

end 


function PeierlsIntegral(A::AbstractVector, B::AbstractVector,
												 ::Val{:y})::Float64
	
	#	Landau Gauge A = [0, x, 0]
	
	(A[1]+B[1])*(B[2]-A[2])

end

function PeierlsPhaseFactor(phi::Real, uc_area::Real, args...
														)::Union{ComplexF64,Function}

	PeierlsPhaseFactor(phi/uc_area, args...)

end


function PeierlsPhaseFactor(B::Real, 
														transl::Union{AbstractString,Char},
														args...
														)::Union{ComplexF64,Function}

	PeierlsPhaseFactor(B, Symbol(lowercase(translation)), args...)

end 



function PeierlsPhaseFactor(B::Real, translation::Symbol,
														)::Union{ComplexF64,Function}

	isapprox(B,0,atol=1e-14) && return 1.0 + 0.0im

	vt = Val(translation)
	
	e = exp(0.5im*pi*B)

	return function p_p_f(ri::AbstractVector,rj::AbstractVector)::ComplexF64

		e^PeierlsIntegral(ri, rj, vt)

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
														Peierls = (0.0,:x),
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

  Zero,One = zero(samplehopp),one(samplehopp) 

  toBasis = ToBasis(using_SC_basis,Zero)

	d0 = LA.checksquare(toBasis(samplehopp))

  same = Utils.fSame(dist_tol)

        # ------------------- needed for the hoppings ----------------- #      


  Append!,Sum,Hoppings = TBmodel.Add_Hopping_Terms(d0,hopp_cutoff)
		# value, condition, (function )
  
  maxdist,update = TBmodel.Update_Maxdist(dist_tol, hopp_cutoff)

  
        # -------------- nearest neighbor ------------------------ #      

  if length(dist)>=1

		Append!(Hoppings, param_H[:Hopping], same(dist[1]) ∘ -, 
						toBasis(One, PeierlsPhaseFactor(param_H[:Peierls]...)))

    maxdist = update(maxdist, param_H[:Hopping], dist[1])

  end  



        # -------------- atom bias (lattice imbalance) ----------- #      

#  Hoppings = Append(Hoppings, param_H[:Lattice_Imbalance],
#		(ri,rj) -> same(ri,rj),
#		(ri,rj) -> 2*ri[sl] - sum(values(slcoord)),
#							)

        # ----------------- local potential --------------- #      

	v0, vf = args_local_potential(; param_H...)

	Append!(Hoppings, v0, same(0), toBasis(One,vf))


        # -------------- superconducting pairing ----------------- #      

  if using_SC_basis 

    val, nmax = param_H[:SC_Gap]

    Append!(Hoppings, 1, true, val)

    nr_uc = max(Int64(nmax), nr_uc)

    maxdist = update(maxdist, 1, nmax)

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

  function cond(ri::AbstractVector,rj::AbstractVector)::Bool 

		(isnothing(indomain) || indomain(ri,rj)) && LA.norm(ri-rj)<maxdist

	end
  
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

