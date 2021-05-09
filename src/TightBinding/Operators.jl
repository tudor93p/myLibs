module Operators

import ..LA, ..SpA 

import ..Utils, ..TBmodel, ..Algebra




#===========================================================================#
#
# Special operators usually used
# 
#---------------------------------------------------------------------------#



function Velocity(BlochH, dir=1; dk=pi/100)

	isnothing(dir) | isnothing(BlochH) && return (P;k) -> zeros(size(P,2),1)

	function v(k) 
		
		k1, k2 = copy(k), copy(k)

		k1[dir] -= dk/2

		k2[dir] += dk/2

		return (BlochH(k2) - BlochH(k1))/dk

	end

#	(BlochH(k .+ dk*(axes(k,1).==dir)) - BlochH(k))/dk

  return (P;k,kwargs...) -> reshape(LA.diag(P'*v(k)*P),:,1)

end


	# ----- probability of localization on individual atoms ----- #

function LDOS(Op=[1];kwargs...)

  Operator(Op,sum_atoms=false;kwargs...)


end


	# ----- inverse participation ratio ----- #


function IPR(;kwargs...)

  ldos = LDOS(;kwargs...)

  return (P;kwargs...) -> 1 ./ sum(ldos(P).^2,dims=2)
		# special function, cannot be constructed from Operator(...)
		# = 1 if psi localized on a single site, 
		# = nr_atoms for psi spread on all atoms
end



	# ----- expectation of the position operator ----- #


function Decode_PositionExpVal(str;err=false)
	
	good(A...) = any(isequal(str), A)
									 
	for (d,A) in enumerate(["X","Y","Z"])
	
		good(A) && return d, identity
	
		good("|$A|") && return d, abs
	
		good(A*"2", "|$A|2") && return d, abs2 
		
		good(A*"3") && return d, x->x^3

		good("|$A|3") && return d, x->abs(x)^3

		good(A*"4", "|$A|4") && return d, x->x^4 
	
	end 

	err && error("'$str' not understood")

	return nothing,nothing

end 



function Partial_PositExpect_fromLDOS(
					Prob::AbstractMatrix, atoms::AbstractMatrix, str::AbstractString;
					kwargs...)

	Partial_PositExpect_fromLDOS(Prob, atoms,
															 Decode_PositionExpVal(str, err=true)...;
															 kwargs...)

end







function Partial_PositExpect_fromLDOS(
								Prob::AbstractMatrix, atoms::AbstractMatrix,
								dim::Int, f::Function=identity; 
								convolute=false, delta=0.3, vsdim=2,
								kwargs...)

# P[energy, atom] = LDOS at energy, on atom 
# or: P[wf_index, atom] = localiz.prob on atom of a certain wf

#	A = the coordinates 'dim'. B =  the perpendicular



	dim_perp = ifelse(dim==1, 2, ifelse(dim==2, 1, nothing)) 


	function out(get_ylim_Pf, param_B)

		ylim, Prob_and_fA = get_ylim_Pf()

		exp_vals = zeros(length(param_B), size(Prob,1)) 
	
		for i in axes(exp_vals,1)

			Pi,fi = Prob_and_fA(i)

			@warn "WFs on columns" 

			exp_vals[i,:] = Algebra.Normalize(Pi,1; dim=vsdim)*fi

		end

		#	out[energy,some_y] = 	(sum over i in atoms_at_y)(P[energy, i] * f(x[i]))
		#												/(sum ...) P[energy, i]
	

#		exp_vals = Algebra.Normalize_Columns(exp_vals, 1)



		return (param_B, exp_vals, ["x","y","z"][dim_perp], extrema(ylim))

	end


	sel(X,i) = selectdim(X, [2,1][vsdim], i)

	if !convolute


		uniq_B, inds_B = Utils.Unique(sel(atoms, dim_perp),
																	inds="all", sorted=true)

		return out(uniq_B) do 

			f_A = f.(sel(atoms, dim))

			return (f_A, function(i) 

						inds_B[i] |> I -> (sel(Prob, I), f_A[I])

									end)

		end


	end  



	


	(A,dense_A,w_A), (B,dense_B,w_B) = map((dim,dim_perp)) do d
	
		v = sel(atoms, d)
	
		u = Utils.Unique(v, sorted=true)

		m,M = extrema(u)

#		m,M = [m,M] .+ (d==dim)*[1,-1]*1/4*(M-m) # disregard edges
	
		return (v, range(m, M, length=150), minimum(diff(u))*delta)
	
	end


	return out(dense_B) do 

		convProb = Prob*Algebra.get_CombinedDistribution(
												(A, B),
												(repeat(dense_A, outer=length(dense_B)),
																 repeat(dense_B, inner=length(dense_A))),
												(w_A, w_B);
												weights="Lorentzian", normalize=false)

		LI = LinearIndices((axes(dense_A,vsdim), axes(dense_B,vsdim)))

		f_A = f.(dense_A)

		return (f_A, function(i)
					s = LI[CartesianIndices((axes(dense_A,vsd), i))][:]

					(sel(convProb,s), f_A)

				end)
	end

#	return out(get_ylim_Pf, dense_B)




#	L = Lorentzian with withd dB/5 or dA/5
#	P(R) = mapreduce(+, axes(atoms,1)) do i 
#				mapreduce(a->L(a[1]-a[2]), *, zip(R,atoms[i,:]), init=LDOS[i])
#			end

			
end





function Position_Expectation(axis::Int64, Rs::AbstractMatrix{Float64};
															fpos::Function=identity,dim=2,
															nr_at=size(Rs,dim), kwargs...)


  Operator(map(fpos, selectdim(Rs, [2,1][dim], axis)),
					 diag="orbitals", nr_at=nr_at; kwargs...)

end



	# ----- local charge operator ----- #


function Charge(Op_UC,atom;kwargs...)

  Operator(Op_UC,nr_orb=size(Op_UC,1),purpose="matrixelement",acts_on_atoms=vcat(atom);kwargs...)

end



function Norm(;kwargs...)

	Operator([1];kwargs...)

end


function Overlap(;kwargs...)

  Operator([1],purpose="matrixelement";kwargs...)

end


#===========================================================================#
#
# Process input and construct the operator with a specific structure
# 
#---------------------------------------------------------------------------#


function Operator(Op::AbstractArray; purpose="expectation", kwargs...)

  Op,(nr_at,nr_orb,size_H),diag = Understand_OperatorInput(Op;kwargs...)
		# checks if the size of Op coincides with 
		# 	the number of orbitals or atoms 
		#	or the Hamiltonian size 
		# 	and the consistency with 'diag'
		#	diag can be ["orbitals","atoms","all","no"]


	p = findfirst(isequal(purpose), ["expectation","matrixelement"]) 
	
	isnothing(p) && error("'purpose' should be either 'expectation' or 'matrixelement', not '$purpose'")



  function choose(operators)

		out = operators[p]
	
		!isnothing(out) && return out 

		error("This operator cannot be used for '$purpose'")

  end



#  println(diag,"  ",purpose,"  ",nr_orb,Op)


  if diag == "no" 
 
#    println("full") 
#    exit()
    return Operator_Full(Op) |> choose

  elseif diag =="all"
 
#    println("1x1")
#    exit()
    return Operator_aNumber(Op, nr_orb=nr_orb, nr_at=nr_at; kwargs...) |> choose
 
 
  elseif diag in ["orbitals", "atoms"]

#    println("decomposable")  
#    exit()

    return Operator_Decomposable(Op, diag=diag, nr_orb=nr_orb, nr_at=nr_at; kwargs...) |> choose


  end


end


#===========================================================================#
#
# Operator for the cases 	O = O_{atoms}⊗ 1_{orbitals} 
#			or 	O = 1_{atoms}⊗ O_{orbitals}	
# 
#---------------------------------------------------------------------------#

function Operator_Full(Op::AbstractArray)


  return [(P;kwargs...) -> ExpectationValues_FullOp(P,Op),
	(Pl,Pr=Pl;kwargs...) -> MatrixElements_FullOp(Pl,Pr,Op)]


end


function Operator_Decomposable(Op::AbstractArray; 
		diag=nothing,
		sum_atoms::Bool=true, sum_orbitals::Bool=true,
		acts_on=nothing,
		nr_at=nothing, nr_orb=nothing, Kwargs...)

  if isnothing(diag) || !in(diag,["orbitals","atoms"])
    error("Please specify if the operator is diagonal in atoms or orbitals")
  end

  


  function get_acts_on(acts_on_, diag_, nr_at_, nr_orb_)
  
    !isnothing(acts_on_) && return acts_on_
  
    diag_=="atoms" && !isnothing(nr_at_) && return 1:nr_at_
  
    diag_=="orbitals" && !isnothing(nr_orb_) && return 1:nr_orb_
  
    error("It's unclear on what the operator acts")
  
  end


  function argsinds(acts_on_,diag_,nr_at_,nr_orb_)

    diag_=="atoms" && return (1:nr_orb_, acts_on_)

    diag_=="orbitals" && return (acts_on_, 1:nr_at_)

  end


  function repeatf(acts_on_, diag_, nr_at_, nr_orb_)

    if diag_=="atoms" 

      return o -> Repeat_Operator_ManyAtoms(o, length(acts_on_), nr_orb_)
     
    elseif diag_=="orbitals" 

      return o -> Repeat_Operator_ManyOrbitals(o,nr_at_,length(acts_on_))

    end
  end




  acts_on = get_acts_on(acts_on, diag, nr_at, nr_orb)
		# if acts_on not specified, all orbitals/atoms are assumed



  inds = TBmodel.Hamilt_indices_all(argsinds(acts_on,diag,nr_at,nr_orb)..., nr_orb; iter=diag)
	# for each item in iter: gives the indices of wf components
	#		on which the operator acts

  all_inds = sort(vcat(inds...))
		# indices for all the wf components on which the operator acts


  Repeat_Many = repeatf(acts_on, diag, nr_at, nr_orb)
		# repeat the operator for several atoms/orbitals 



  if (diag=="atoms" && sum_atoms) | (diag=="orbitals" && sum_orbitals)
#  if sum_action

    FullOp = ndims(Op)==1 ? LA.diag(Repeat_Many(LA.diagm(0=>Op))) : Repeat_Many(Op)




    return [(P;kwargs...) -> ExpectationValues(P, FullOp, all_inds), 

	(Pl,Pr=Pl;kwargs...) -> MatrixElements(Pl, Pr, FullOp, all_inds)]

  end


  return [(P;kwargs...) -> ExpectationValues(P,Op,inds,iterate=true),
	[(Pl,Pr=Pl;kwargs...) -> MatrixElements_DiagonalOp(Pl,Pr,Op,i) for i in inds]
			]

end



#===========================================================================#
#
# Operator for the case O = 1_{atoms}⊗ 1_{orbitals}
# 
#---------------------------------------------------------------------------#

function Operator_aNumber(Op::AbstractArray;
		sum_atoms::Bool=true,acts_on_atoms=nothing,
		sum_orbitals::Bool=true,acts_on_orbs=nothing,
		nr_at=nothing,nr_orb=nothing,Kwargs...)



  if isnothing(acts_on_orbs) && !isnothing(nr_orb)
    acts_on_orbs = 1:nr_orb
  end

  if isnothing(acts_on_atoms) && !isnothing(nr_at)
    acts_on_atoms = 1:nr_at
  end


  
  acts_on,iter_name = [acts_on_atoms,acts_on_orbs],["atoms","orbitals"]

  good_iter = (!isnothing).(acts_on) 

  acts_on,iter_name = acts_on[good_iter],iter_name[good_iter]

  shorter_iterator = iter_name[argmin(length.(acts_on))]

  inds = TBmodel.Hamilt_indices_all(acts_on_orbs,acts_on_atoms,nr_orb;iter=shorter_iterator)

  all_inds = sort(vcat(inds...))

  sum_iterator = ["atoms","orbitals"][[sum_atoms,sum_orbitals]]




	# if no sum must be performed
  if length(sum_iterator)==0
     return [(P;kwargs...) -> Op[1]*transpose(abs2.(P[all_inds,:])),
		nothing]
		# matrix element doesn't make sense

	# if both sums must be performed   
  elseif length(sum_iterator)==2

    return [(P;kwargs...) -> Op[1]*transpose(sum(abs2.(P[all_inds,:]),dims=1)),
		(Pl,Pr=Pl;kwargs...)-> Op[1]*Pl[all_inds,:]'*Pr[all_inds,:]
		]

  end




  if shorter_iterator == sum_iterator[1]
    
    return [(P;kwargs...) -> 
		Op[1]*transpose(mapreduce(i->abs2.(P[i,:]),+,inds)),
	nothing] 
  
  else         
    return [(P,kwargs...) -> 
		Op[1]*hcat(map(i->sum(abs2.(P[i,:]),dims=1)[:],inds)...),
	nothing]

  end
end



#===========================================================================#
#
# Various ways to act with operators
# 
#---------------------------------------------------------------------------#

function ExpectationValues(P::AbstractMatrix,M::AbstractVector,inds=Colon(); 
													 iterate=false)

  ExpectationValues_DiagOp(P, M, inds; iterate=iterate) 

end 


function ExpectationValues(P::AbstractMatrix, M::AbstractMatrix, inds=Colon();
													 iterate=false)

	ExpectationValues_FullOp(P, M, inds; iterate=iterate)

end



function MatrixElements(Pl::AbstractMatrix, Pr::AbstractMatrix,
												M::AbstractVector, inds=Colon())

	MatrixElements_DiagOp(Pl,Pr,M,inds)
end 


function MatrixElements(Pl::AbstractMatrix, Pr::AbstractMatrix,
												M::AbstractMatrix, inds=Colon())

	MatrixElements_FullOp(Pl, Pr, M, inds)

end



# ---- 


function ExpectationValues_DiagOp(P::AbstractMatrix, M::AbstractVector,
																	inds=:; dim=2, iterate=false) 

	if !iterate 
		A = abs2.(selectdim(P, [2,1][dim], inds))

		dim==2 && return reshape(M,1,:)*A 
		dim==1 && return A*reshape(M,:,1) 

	end 


	return mapreduce(hcat, inds) do i 
		
		 ExpectationValues_DiagOp(P, M, i)[:] 

	end 

end


function MatrixElements_DiagOp(Pl::AbstractMatrix, Pr::AbstractMatrix,
															 M::AbstractVector, inds=:; dim=2)
	@warn "vectors on columns might not work "

	dim==2 && return Pl[:,inds]'*(reshape(M,1,:).*Pr[:,inds])


	dim ==1 && return Pl[inds,:]'*(reshape(M,:,1).*Pr[inds,:])

end

function MatrixElements_FullOp(Pl::AbstractMatrix,Pr::AbstractMatrix,M::AbstractMatrix,inds=:; dim=2)

#	selectdim(Pl, dim, inds)

#	Utils.CombsOfVecs(M,selectdim(Pr, dim, inds))


  dim==1 && return Pl[inds,:]'*M*Pr[inds,:]

	dim==2 && return Pl[:,inds]'*M*Pr[:,inds]

end




function ExpectationValues_FullOp(P::AbstractMatrix, M::AbstractMatrix,
																	inds=:; iterate=false)
@warn "vectors on columns might now work"

  !iterate && return reshape(LA.diag(MatrixElements_FullOp(P,P,M,inds)),:,1)

  return hcat(map(i->ExpectationValues_FullOp(P,M,i)[:],inds)...)


# diag(P'Op*P) is faster, although <psi_i|Op|psi_j> are computed as well, with i!=j
end





#===========================================================================#
#
# Repeat an orbital-operator for several atoms
#	or an atom-operator for several orbitals
#
#---------------------------------------------------------------------------#


function Repeat_Operator_ManyAtoms(Op,nr_at,nr_orb)
	# Op simulates a Hamiltonian with nr_orb orbitals and one atom
	# 	will be reproduced for nr_at atoms

  i,j,v = SpA.findnz(SpA.sparse(Op))

  indsvals(a) = hcat(TBmodel.Hamilt_indices.([i,j],a,nr_orb)...,v)

  return TBmodel.indsvals_to_Tm(vcat(indsvals.(1:nr_at)...),nr_at*nr_orb)

end


function Repeat_Operator_ManyOrbitals(Op,nr_at,nr_orb)
	# Op simulates a Hamiltonian with nr_at atoms and one orbital
	# 	will be reproduced for nr_orb orbitals




  ai,aj,v = SpA.findnz(SpA.sparse(Op))


  indsvals(o) = hcat(TBmodel.Hamilt_indices.(o,[ai,aj],nr_orb)...,v)


  return TBmodel.indsvals_to_Tm(vcat(indsvals.(1:nr_orb)...),nr_at*nr_orb)

end

#===========================================================================#
#
# Understand how the operator can be simplified and provide the information
# 
#---------------------------------------------------------------------------#



function Understand_OperatorInput(A::AbstractVecOrMat; 
																	diag=nothing, kwargs...) 
 
  if ndims(A)==2 

    if size(A,1) != size(A,2)

      A = reshape(A,:)

		elseif isapprox(A, LA.diagm(0=>LA.diag(A)), atol=1e-6)

      A = LA.diag(A)

    end
  
  end


  dim = size(A,1)

  numbers, enough_data = get_nratorbH(;kwargs...)



  conditions = [prod(size(A))==1;
      		(n->!isnothing(n) && n==dim).(numbers)...]


  if any(conditions)

    meaning = ["all","orbitals","atoms","no"][argmax(conditions)]

    if meaning != Utils.Assign_Value(diag,meaning)

      error("Strange operator input.")

    end
   
    return A,numbers,meaning

  end

  return A,numbers,diag
 



#  conditions = (n->!isnothing(n) && n==dim).(numbers)
#  conditions = [conditions;[all(isnothing.(numbers)),prod(size(A))==1]]
#  meanings = ["orbitals","atoms","no","no","all"][conditions]


#  length(meanings) == 0 && error("Strange operator input.")


#  (isnothing(diag) || diag==meanings[1]) && return A,numbers,meanings[1]


end


function get_nratorbH(nr_at::Integer, nr_orb::Integer, size_H::Integer)

	nr_at*nr_orb==size_H && return (nr_at, nr_orb, size_H), true

	error("Wrong number of atoms/orbitals")

end


function get_nratorbH(nr_at::Nothing, nr_orb::Integer, size_H::Integer)

	get_nratorbH(div(size_H,nr_orb), nr_orb, size_H)

end

function get_nratorbH(nr_at::Integer, nr_orb::Nothing, size_H::Integer)

	get_nratorbH(nr_at, div(size_H,nr_at), size_H)

end

function get_nratorbH(nr_at::Integer, nr_orb::Integer, size_H::Nothing)

	get_nratorbH(nr_at, nr_orb, nr_at*nr_orb)

end


function get_nratorbH(;nr_at=nothing,nr_orb=nothing,size_H=nothing,kwargs...)

	(nr_at,nr_orb,size_H) |> function (args)
		
		count(!isnothing, args)<2 && return args, false

		return get_nratorbH(args...)

	end

end



#===========================================================================#
#
# Normalize wavefunctions on columns
# 
#---------------------------------------------------------------------------#





#===========================================================================#
#
# Trace over orbitals/atoms with some projection operator
# 
#---------------------------------------------------------------------------#

get_diag(A::AbstractMatrix) = LA.diag(A)
get_diag(A::AbstractVector) = A


function Trace(what, Op=[1]; sum_up=false, kwargs...)

#	@warn "Using Trace in Operators with:" Op kwargs

	oa = ["orbitals","atoms"]

	what in oa || error("The first argument must be either of "*string(oa))


	S = size(Op,1)

	if sum_up && (S==1 || all(isapprox.(Op[1], Op[2:end], atol=1e-8))) 
		
		return A -> Op[1]*((typeof(A)<:AbstractMatrix) ? LA.tr : sum)(A)

	end


	(nr_at, nr_orb, size_H), enough_data = get_nratorbH(; kwargs...)


	!enough_data && return function (A)

								Trace(what, Op; sum_up=sum_up, size_H=size(A,1), kwargs...)(A)

									end

	Op = get_diag(Op)


	function check(A) 
		
		@assert size(A,1)==size_H 

		return A

	end


#	nr_at==1 && what=="atoms"





	S==size_H && sum_up && return A -> sum(Op.*get_diag(check(A)))


	what2 = oa[findfirst(!isequal(what), oa)]


	inds = TBmodel.Hamilt_indices_all(1:nr_orb,1:nr_at; iter=what2)


	iter(M) = get_diag(M) |> m -> (m[i] for i in inds)


	S==1 && !sum_up && return A -> Op[1]* sum.(iter(check(A)))


	S==size_H && return A -> [sum(oi.*ai) for (oi,ai) in zip(iter(Op),iter(check(A)))]





	sum_up_ = sum_up ? sum : identity


	nr = Dict("orbitals"=>nr_orb, "atoms"=>nr_at)



	S==nr[what] && return A -> sum_up_([sum(Op.*ai) for ai in iter(check(A))])

	S==nr[what2] && return A -> sum_up_(Op.*sum.(iter(check(A))))


	error("Something went wrong")



end










#############################################################################









end
