module TBmodel
#############################################################################


import ..LA, ..SpA

import ..Utils, ..Algebra, ..Lattices


#===========================================================================#
#
# Calculates 	--> the Hamiltonian of the layer parallel to a given surface
#							--> the perpendicular hopping Hamiltonian
#
#---------------------------------------------------------------------------#


function default_surface(Latt::Lattices.Lattice)::Int 
	
	@assert Lattices.LattDim(Latt)>0

	return 1

end 



function args_BHPP(Latt::Lattices.Lattice, 
									 arg1_BH::Function,
									 surface::Int=default_surface(Latt),
									 )::Tuple{Lattices.Lattice, Function, Dict, Int}

	(Latt, arg1_BH, Dict(), surface) 

end 


function args_BHPP(Latt::Lattices.Lattice, arg2_BH::AbstractDict, args...
									 )::Tuple{Lattices.Lattice, Function, Dict, Int}

	args_BHPP(Latt, Lattices.NearbyUCs, arg2_BH, args...)

end 

function args_BHPP(Latt::Lattices.Lattice,
									 surface::Int=default_surface(Latt),
									 )::Tuple{Lattices.Lattice, Function, Dict, Int}

	(Latt, Lattices.NearbyUCs, Dict(), surface)

end 


function args_BHPP(Latt::Lattices.Lattice,
									 arg1_BH::Function,
									 arg2_BH::AbstractDict,
									 surface::Int=default_surface(Latt),
									 )::Tuple{Lattices.Lattice, Function, Dict, Int}

	(Latt, arg1_BH, arg2_BH, surface)

end 









function BlochHamilt_Parallel_(Latt::Lattices.Lattice,
															arg1_BH::Function,
															arg2_BH::AbstractDict,
															surface::Int;
															argH="phi"
															)::Function
 
	# Remove the specified dimension "surface" => dimension = d-1
  
	Bloch_Hamilt(arg1_BH(Lattices.ReduceDim(Latt, surface));
							 arg2_BH..., argH=argH)

	# Can be a function of  ---- k_parallel, a vector of length d-1
	# 											---- the Bloch phase phi = k_parallel dot a
	#
  
end

function BlochHamilt_Parallel(args...; kwargs...)::Function 

	BlochHamilt_Parallel_(args_BHPP(args...)...; kwargs...)

end 









#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#





function iter_matrix(A::AbstractMatrix{T})::Vector{Vector{T}} where T<:Number

	[collect(a[:]) for a in eachcol(A)]

end 


function iter_matrix(A::AbstractVector{T})::Vector{T} where T<:Number 

	collect(A)

end 

	
function iter_matrix(A::AbstractVector{<:AbstractVector{<:T}})::Vector{Vector{T}} where T<:Number

	[collect(a[:]) for a in A]

end 




function BlochHamilt_Perpendicular(args...)::Dict 
	
	BlochHamilt_Perpendicular_(args_BHPP(args...)...)
	
end 

 
function BlochHamilt_Perpendicular_(Latt::Lattices.Lattice,
																		arg1_BH::Function, 
																		arg2_BH::AbstractDict,
																		surface::Int,
																		)::Dict

	
        # Remove from Latt all dimensions != "surface" 
				# the resulting lattice has dimension 1
	
	TBL = arg1_BH(Lattices.KeepDim(Latt, surface))
	
	ms = [length(m)==1 ? first(m) : error() for m in iter_matrix(TBL[1])]


	Rms = Utils.VecAsMat.(iter_matrix(TBL[2]), 2)


	AtomsUC = TBL[3] 


	for sgn in [+1,-1]

		useful_ms = ms*sgn .> 0 

		any(useful_ms) || continue 
		
		return Dict(map(sgn * ms[useful_ms], sgn * Rms[useful_ms]) do m,R

			m => HoppingMatrix(AtomsUC, AtomsUC .+ R; arg2_BH...)

		end)

	end

	return Dict()

end

function BlochHamilt_ParallPerp(args...)::Tuple{Function, Dict}
	
	Utils.invmap(args_BHPP(args...), (BlochHamilt_Parallel_,
																		BlochHamilt_Perpendicular_))

end 




#===========================================================================#
#
#		Defines the arrangement of the Hamiltonian:
#
#					[ all orbitals of atom 1, all orbitals of atom 2, ... ]
#				
#---------------------------------------------------------------------------#

function Hamilt_indices(orbital::Int, atom::Int, nr_orbitals::Int)::Int

	orbital + (atom - 1) * nr_orbitals

end

function Hamilt_indices(orbitals::AbstractVector{<:Int}, 
												atom::Int, nr_orbitals::Int)::Vector{Int}

	Hamilt_indices.(orbitals, atom, nr_orbitals)

end



function Hamilt_indices(orbital::Int,
												atoms::AbstractVector{<:Int}, 
												nr_orbitals::Int)::Vector{Int}

	Hamilt_indices.(orbital, atoms, nr_orbitals)

end



function Hamilt_indices_all(orbs::Utils.List, atoms::Utils.List,
														d0::Int=length(orbs);
														iter::AbstractString="atoms",flat::Bool=false
														)::Vector

	iter=="orbitals" && flat==true && @warn "The function returns flat  Hamiltonian indices by iterating over orbitals. In multi-orbital system, the order will *not* be consistent with the way the Hamiltonian is built and unexpected results will occur.\n"

#  f(x) = flat ? vcat(x...) : x

#  iter == "orbitals"  && return f(map(i->Hamilt_indices(i,atoms,d0),orbs))
  
	inds_generator = if iter=="orbitals"  
		
		(Hamilt_indices(i,atoms,d0) for i in orbs)

										elseif iter=="atoms"

		(Hamilt_indices(orbs,i,d0) for i in atoms) 

										end 


	return flat ? vcat(inds_generator...) : collect(inds_generator)


  
#  iter == "atoms" && return f(map(i->Hamilt_indices(orbs,i,d0),atoms))

  error("'iter' should be either 'orbitals' or 'atoms'")

end


#function combineOperators_AtOrb(op_at::SpA.AbstractSparseMatrix{Ta},
#																op_orb::AbstractMatrix{To},
##											 )#::SpA.SparseMatrixCSC{promote_type(Ta,To)
#												#									} 
#												) where {Ta<:Number, To<:Number}
#	println("sparse at")
#
#	nr_at = LA.checksquare(op_at)
#	
#	nr_orb = LA.checksquare(op_orb)
#
#	size_H = nr_at*nr_orb
#
#	FullOp = SpA.spzeros(promote_type(Ta,To), size_H, size_H)
#
#
#	inds = Hamilt_indices_all(1:nr_orb, 1:nr_at; iter="atoms")
#
#
#	for (i, j, v_ij) in zip(SpA.findnz(op_at)...)
#
#		FullOp[inds[i], inds[j]] .= op_orb * v_ij 
#
#	end 
#
#	return FullOp
#
#end   
#
#
#
#function combineOperators_AtOrb(op_at::AbstractMatrix{Ta},
#																op_orb::SpA.AbstractSparseMatrix{To},
##											 )
#	#::SpA.SparseMatrixCSC{promote_type(Ta,To)
#												#									} 
#												) where {Ta<:Number, To<:Number}
#
#
#	println("sparse orb")
#
#	nr_at = LA.checksquare(op_at)
#	
#	nr_orb = LA.checksquare(op_orb)
#
#	size_H = nr_at*nr_orb
#
#	FullOp = SpA.spzeros(promote_type(Ta,To), size_H, size_H)
#
#
#	inds = Hamilt_indices_all(1:nr_orb, 1:nr_at; iter="orbitals")
#
##	I,J,V = SpA.findnz(op_orb)
#
#	for (i, j, v_ij) in zip(SpA.findnz(op_orb)...)
#
#		FullOp[inds[i], inds[j]] .= op_at * v_ij 
#
#	end 
#
#	return FullOp
#
##	return kron!(FullOp, op_at, op_orb)
#
#end   

function combineOperators_AtOrb(op_at::AbstractMatrix{<:Number},
																op_orb::AbstractMatrix{<:Number},
																)::AbstractMatrix{<:Number} #%where {Ta<:Number, To<:Number}

	kron(op_at, op_orb)

#	size_H = LA.checksquare(op_at)*LA.checksquare(op_orb)
#
#	FullOp = zeros(promote_type(Ta,To), size_H, size_H)
#
#	return kron!(FullOp, op_at, op_orb)

end  



#===========================================================================#
#
#		Given the hopping matrix v between atoms i and j,
#
#				--> if v is a (1x1) matrix : return (i,j,v)
#				--> if v is a (nxn) matrix : return (i1.j1,v[1]), (i2,j2,v[2]), ...,
#								i1,j1 etc. are the indices of the final Hamiltonian
#
#---------------------------------------------------------------------------#

function flat_indsvals(i::Int64, j::Int64, v::Number, d0::Int64=1)::Matrix

	flat_indsvals(i, j, hcat(v), d0)

end 


function flat_indsvals(i::Int64, j::Int64, v::AbstractMatrix, d0::Int64)::Matrix


  I,J,V = SpA.findnz(SpA.sparse(v))

  return hcat(Hamilt_indices.([I,J],[i,j],d0)...,V)

end 



#===========================================================================#
#
#		Transforms the big array [ I | J | V ] to a matrix T, T[I[a],J[a]] = V[a]
#
#---------------------------------------------------------------------------#

function indsvals_to_Tm(indsvals::AbstractMatrix, n::Int64, m::Int64=n
											 )::SpA.SparseMatrixCSC

	types = (Int64, Int64, ComplexF64)

	hopps = (convert(Vector{t}, c) for (t,c) in zip(types, eachcol(indsvals)))

  return SpA.sparse(hopps..., n, m)

end



#===========================================================================#
#
#	Get all nonzero hoppings between 	- atom_i at position Ri 
#																		- and all atoms j at positions Rsj
#
#---------------------------------------------------------------------------#

function get_hopps(atom_i::Int, Ri::AbstractVector, 
									 Rsj::AbstractMatrix, 
									 hopping::Function, hopp_cutoff::Float64, 
									 d0::Int)::Matrix{ComplexF64}

  hopps = hopping.([Ri],eachcol(Rsj))
		# compute the hopping amplitudes from Ri to all Rj-s

	@assert d0==(isa(rand(hopps),Number) ? 1 : LA.checksquare(rand(hopps)))

  atoms_js = findall(LA.norm.(hopps) .> hopp_cutoff)
		# check which amplitudes are nonzero


  isempty(atoms_js) && return zeros(0, 3)

  return vcat(flat_indsvals.([atom_i], atoms_js, hopps[atoms_js], [d0])...)


end


#===========================================================================#
#
# Given the positions of two sets of atoms (Rsi and Rsj)
#		and a hopping function hopping(ri,rj), 
#		it outputs the hopping matrix Tij
#
#---------------------------------------------------------------------------#


function HoppingMatrix(Rsi::AbstractMatrix, Rsj::AbstractMatrix=Rsi;
											 Hopping::Function, 
											 nr_orb::Int=1, 
											 parallel::Bool=false, 
											 hopp_cutoff::Float64=1e-6, 
											 kwargs...
											 )::SpA.SparseMatrixCSC


  indsvals_to_Tm(vcat(

      (parallel ? pmap : map)(enumerate(eachcol(Rsi))) do (i,Ri)
      	# for each atom i at position Ri

				get_hopps(i, Ri, Rsj, Hopping, hopp_cutoff, nr_orb)

	  	# find the atoms where it hopps and the amplitudes

    end...), size(Rsi,2)*nr_orb, size(Rsj,2)*nr_orb)

end

#===========================================================================#
#
# 	Computes all tight-binding matrices
#				for a certain lattice and hopping function
#
#---------------------------------------------------------------------------#

function Compute_Hopping_Matrices(
						Lattice::NTuple{3,AbstractArray};
						dist_tol::Float64=1e-5, parallel::Bool=false,	kwargs...
						)::Tuple{SpA.SparseMatrixCSC, Any}


#  ms, Rms, AtomsUC = Lattice 


  # 'AtomsUC' 	-> the positions of the atoms in the central unit cell
  # 'Rms'	-> the positions of the unit cells around the central one
  # 'ms'	-> the integer vectors ????  designating those unit cells
  	#	-> one ???? vector on each column of the matrix ms
	#	-> ???? has the same length as the number of lattice vectors
	# 	-> each UC in the list is ????_???? = ??_j (m_j ??? ????_j)

  # kwargs include
  # Hopping 	-> the function f(ri,rj) which gives the hopping amplitude
  #	 		where the atoms in the central cell can hop;
  # hopp_cutoff -> the cutoff for the hopping amplitudes
  # nr_orb -> the size of the matrix f(ri,rj) 
  # ...

	AtomsUC = Lattice[3] 

	ms,Rms = map(iter_matrix, Lattice[1:2])


	# if few procs, Rms will be distributed. Otherwise -> atoms

 	Tms = (parallel && nprocs()<=length(Rms) ? pmap : map)(Rms) do Rm

    HoppingMatrix(AtomsUC,	Rm .+ AtomsUC;
													kwargs...,
													parallel = parallel && nprocs()>length(Rms)
      										)
  end


  nonzeroT = findall(SpA.nnz.(Tms) .> 0)
		# check which T_????  contain only zeros


	intra = findfirst(m -> m==-m, ms)
		# corresponds to the intra-cell hopping

  inter = nonzeroT[nonzeroT .!= intra]
		# correspond to the non-zero inter-cell hopping

	isempty(inter) && return Tms[intra], nothing 

	return Tms[intra], Tuple.(getindex.((ms, Rms, Tms), (inter,)))


end


#===========================================================================#
#
# Compute the TB matrices, (save them) and return the Hamiltonian
#
#---------------------------------------------------------------------------#

function Bloch_Hamilt(Lattice::NTuple{3,AbstractArray};
											argH::AbstractString="", 
											savefile::AbstractString="",	
											Hopping...)::Function
			
  data = Compute_Hopping_Matrices(Lattice; Hopping...)

  isempty(savefile) || error("Save_Hamilt_Data not yet rewritten")

	# Save_Hamilt_Data(savefile,mRT...)

  return Assemble_Bloch_Hamiltonian(data...; argH=argH)

end


function Bloch_Velocity(Lattice::NTuple{3,AbstractArray}, dir...;
												argH::AbstractString="",
												Hopping...)::Function
			
  intra, inter = Compute_Hopping_Matrices(Lattice; Hopping...)

	@assert !isnothing(inter) "Local Hamiltonian"

  return Assemble_Velocity_Operator(inter, dir...; argH=argH)

end


function Bloch_FullVelocity(Lattice::NTuple{3,AbstractArray};
														argH::AbstractString="",
														Hopping...)::Vector{Function}
			
	intra, inter = Compute_Hopping_Matrices(Lattice; Hopping...)

	@assert !isnothing(inter) "Local Hamiltonian"
#	isnothing(inter) && return Function[]

#		ms,Rms = map(iter_matrix, Lattice[1:2])
#
#		mRT = ((Int[],), (Rms[1],), (zero(intra),))
#
#		return Assemble_FullVelocity_Operator(mRT, argH=argH)
#
#	end 



  return Assemble_FullVelocity_Operator(inter; argH=argH)

end



function Bloch_NVelocity(Lattice::NTuple{3,AbstractArray}, dir...;
												 kwargs...)::Function

	Bloch_NVelocity(Bloch_Hamilt(Lattice; kwargs...), dir...; kwargs...)

end 

function Bloch_NVelocity(BlochH::Function;
												 dk::Float64=pi*1e-6, kwargs...
												 )::Function

	function v(k::AbstractVector{<:Real})::AbstractMatrix{<:Number}

		@assert length(k)==1

		(BlochH(k[1]+dk/2) - BlochH(k[1]-dk/2))/dk

	end 

end 


function Bloch_NVelocity(BlochH::Function, dir::Int;
												 dk::Float64=pi*1e-6, kwargs...
												 )::Function
	
	@assert 1<=dir<=3

	D = setindex!(zeros(3), dk/2, dir)


	return function v(k::AbstractVector{<:Real})::AbstractMatrix{<:Number}

		@assert 1<=dir<=length(k) 

		(BlochH(k+D[axes(k,1)])-BlochH(k-D[axes(k,1)]))/dk

	end

end


function Bloch_Hamilt(savefile::String="";argH::String="")::Function

  error("Read_Hamilt_Data not yet rewritten")

  return Assemble_Bloch_Hamiltonian(Read_Hamilt_Data(savefile)...,argH=argH)

end


#===========================================================================#
#
# Build a tight-binding Hamiltonian H(k) based on hopping matrices
#
#---------------------------------------------------------------------------#


function Assemble_Bloch_Hamiltonian(intra::AbstractMatrix, 
																		inter::Nothing=nothing; 
																		kwargs...)::Function

	# if this is a finite system or the Hamiltonian is simply local
	
	H_no_k(args...) = intra 

end 

function iter_k_or_phi(ms::NTuple{N, Union{Int,<:AbstractVector{Int}}},
											 Rms::NTuple{N, <:AbstractVector{Float64}};
											 argH::AbstractString="") where N

	@assert !any(m->  m==-m, ms)
	
	[Rms,ms][findfirst(get_argH(argH, ms, Rms).==["k","phi"])]

end 

function BlochPhase(k,R)::ComplexF64

	exp(im * LA.dot(k,R))

end 

function Assemble_Bloch_Hamiltonian(intra::AbstractMatrix, (ms, Rms, Tms);
																		kwargs...)::Function

#		phases = (Algebra.dot(k,R) for R in iter2) 
					# Bloch phase ???^(i ???? ??? ???? ), 


	iter2 = iter_k_or_phi(ms, Rms; kwargs...)

	return function h(k)::AbstractMatrix{<:ComplexF64}

		h1 = mapreduce(+, Tms, iter2) do T,R

			T * BlochPhase(k, R)

		end 

    return intra + h1 + h1'

  end


end


function Assemble_Velocity_Operator(Ts::NTuple{N, <:AbstractMatrix{<:Number}},
																		Rs::NTuple{N, Union{Int, <:AbstractVector{<:Real}}},
																		dir::Int)::Function where N 

	function v_dir(k)::AbstractMatrix{<:ComplexF64}

		v1 = mapreduce(+, Ts, Rs) do T,R 

			T * R[dir] * BlochPhase(k, R) 

		end 

		return im*(v1-v1')

	end 

end 



function Assemble_Velocity_Operator((ms,Rms,Tms), dir::Int=0; 
																		kwargs...)::Function

	iter2 = iter_k_or_phi(ms, Rms; kwargs...)

	dir += (dir==0 && length(iter2[1])==1)

	@assert 1<=dir<=length(iter2[1])

	return Assemble_Velocity_Operator(Tms, iter2, dir)

end 




function Assemble_FullVelocity_Operator((ms,Rms,Tms); kwargs...
																				)::Vector{Function}

	iter2 = iter_k_or_phi(ms, Rms; kwargs...)

	all_dirs = 1:length(iter2[1])

	return [Assemble_Velocity_Operator(Tms, iter2, dir) for dir in all_dirs]

end 





#===========================================================================#
#
#		Helper function for deciding what argument the Hamiltonian gets
#
#---------------------------------------------------------------------------#

function get_argH(argH::AbstractString, 
									ms::NTuple{N, Union{Int, AbstractVector{Int}}},
									Rms::NTuple{N, AbstractVector{Float64}},
									)::String where N 
  

  args = ["k","phi"]

  argH in args && return argH

	@assert isempty(argH) "** Error in Bloch_Hamiltonian: argH=$argH not understood. It should be an empty string or either of $args **"

	# it should be physical 'k' if the length of the lattice vectors,
	#	 is equal to their number, and Bloch phase 'phi' otherwise

	return length(ms[1])==length(Rms[1]) ? args[1] : args[2]

end


#===========================================================================#
#
# Save the tight-binding Hamiltonian
#
#---------------------------------------------------------------------------#

function Save_Hamilt_Data(filename,ms,Rms,Tms)

	error("Outdated function!")

  open(filename,"w") do f

    DelimitedFiles.writedlm(f,[[size(ms,1),size(ms[1],1)],size(Tms[1])])
  
    DelimitedFiles.writedlm(f,ms)

    DelimitedFiles.writedlm(f,Rms)
 
    (M->[DelimitedFiles.writedlm(f,[item]) for item in SpA.findnz(M)]).(Tms)
 

  end
end




function Read_Hamilt_Data(filename)


	error("Outdated function!")

  isfile(filename) || error("Hamiltonian files not found")


  open(filename,"r") do f

    readnrs(typ) = parse.(typ,split(strip(readline(f),'\n'),'\t'))

    (nr_m,dim_m),dim_Tm = readnrs.([Int64,Int64])

    matrix(t) = readnrs.(repeat([t],nr_m))   
 
    Tm(i) = SpA.sparse(readnrs.([Int64,Int64,Complex{Float64}])...,dim_Tm...)


    return (matrix.([Int64,Float64])...,Tm.(1:nr_m))

  end




end



#===========================================================================#
#
#		Keep track of the largest bond (while adding hopping terms)
#
#---------------------------------------------------------------------------#



function Update_Maxdist(dist_tol::T, hopp_cutoff::Real
											 )::Tuple{T,Function} where T<:Real

  function update(maxdist::Real, param::Real, d::Real)::Float64

    abs(param)>hopp_cutoff ? max(maxdist,d)+2dist_tol : maxdist

  end

  return 2dist_tol,update

end

function Update_Maxdist!(dist_tol::T, hopp_cutoff::Real
												)::Tuple{Vector{T},Function} where T<:Real

	function update!(md::AbstractVector{<:Real}, param::Real, d::Real 
									 )::Nothing 

		abs(param)>hopp_cutoff && setindex!(md, max(md[1],d)+2dist_tol, 1)

		return 

  end

	return [2dist_tol],update!

end



function constHoppTerm(cht::Union{AbstractMatrix{<:Number},Number}
											)::Function 
    
	function f(ri::AbstractVector{<:Real}, rj::AbstractVector{<:Real}
						 )::Matrix{ComplexF64}

		copy(matrixval(cht))

	end 

end 


#===========================================================================#
#
# Initialize hopping terms -- different methods
# 		 
#
#---------------------------------------------------------------------------#

function matrixval(val::Number)::Matrix{ComplexF64}
	
	hcat(val)

end 

function matrixval(val::AbstractMatrix)::Matrix{ComplexF64}

	val 

end 


function matrixval(f::T)::T where T<:Function 

	f

end 


function matrixval(v1::Union{Number,AbstractMatrix},
									 v2::Union{Number,AbstractMatrix})::Matrix{ComplexF64}

	matrixval(v1).*matrixval(v2)

end


function matrixval(v::Union{Number,AbstractMatrix},
									 f::Function)::Function 

	V = matrixval(v) 
	
	return F(args...)::Matrix{ComplexF64} = V .* f(args...) 

end









function Add_Hopping_Terms(d0::Int, hopp_cutoff::Float64=1e-8)

	small = Utils.fSame(hopp_cutoff)(0)

  Zero = zeros(Complex{Float64}, d0, d0)


  #function Append(hoppings::T, value, args...
	#								)::T where T<:AbstractVector{Function}

  #  small(value) && return hoppings

  #  return [hoppings; Hopping_Term(matrixval(value), Zero, args...)]
  #end

  function Append!(hoppings::T,
									 value,args...
									 )::T where T<:AbstractVector{Function}

    small(value) && return hoppings

		return push!(hoppings, Hopping_Term(matrixval(value), Zero, args...))

  end



	function Sum_(hoppings::AbstractVector{Function})::Function

		length(hoppings)==1 && return only(hoppings)

		return function total_hopp(ri::AbstractVector,
															 rj::AbstractVector)::AbstractMatrix 
			
			out = copy(Zero)

			for h in hoppings 

				out += h(ri,rj)

			end 

			return out 

		end 

	end

  function Sum(hoppings::AbstractVector{Function}, condition)::Function

		isempty(hoppings) && return constHoppTerm(Zero)

		return Hopping_Term(matrixval(1), Zero, condition, Sum_(hoppings))

  end



  return Append!, Sum, Function[]


end


	# ------------ General case, matrix ------------ #

function Hopping_Term(value::AbstractMatrix,
											Zero::AbstractMatrix,
											condition::Function,
											fun::Function)::Function


  function f(ri::AbstractVector, rj::AbstractVector)::AbstractMatrix

    condition(ri,rj) ? value .* fun(ri,rj) : Zero

  end

end

	# ------------ Constant matrix if condition  ------------ #

function Hopping_Term(value::AbstractMatrix,
											Zero::AbstractMatrix,
											condition::Function,
											fun::AbstractMatrix)::Function

	function f(ri::AbstractVector,rj::AbstractVector)::AbstractMatrix
	
	  condition(ri,rj) ? value .* fun : Zero
	
	end

end


	# ------------ Always return the hopping, function ------------ #

function Hopping_Term(value::AbstractMatrix,
											Zero::AbstractMatrix,
											condition::Bool,
											fun::Function)::Function

	condition || return constHoppTerm(Zero)

  return function f(ri::AbstractVector, rj::AbstractVector)::AbstractMatrix
  
    value .* fun(ri,rj)

  end

end

	# ------------ Always return the hopping, const ------------ #

function Hopping_Term(value::AbstractMatrix,
											Zero::AbstractMatrix,
											condition::Bool,
											fun::AbstractMatrix)::Function

	condition || return constHoppTerm(Zero)

  return function f(ri::AbstractVector, rj::AbstractVector)::AbstractMatrix
  
    value .* fun

  end

end

	# ------------ Constant hopping, matrix  ------------- #


function Hopping_Term(value::AbstractMatrix,
											Zero::AbstractMatrix,
											condition::Function)::Function

	function f(ri::AbstractVector,rj::AbstractVector)::AbstractMatrix

    condition(ri,rj) ? value : Zero

  end

end



	# ------------ Hopping term from TB file ------------- #

function Hopping_Term_fromTBfile(filename)

  #TBmodel.Save_Hamilt_Data(filename,ms,Rms,Tms) #was done before

  Hopping_Term_fromTBmatr(Read_Hamilt_Data(filename)...)

end


function Hopping_Term_fromTBmatr(ms,Rs,Ts)::Function

  function condition_value_maxdist(same,ax=axes(Rs,1))


    gooddist(ri,rj) = same.(eachcol(Rs[ax,:]),[(ri .- rj)[ax]])

    cond(ri,rj) = any(gooddist(ri,rj))

    val(ri,rj)  = Ts[argmax(gooddist(ri,rj))]

    maxdist = maximum(LA.norm, eachcol(Rs[ax,:]))*1.05
   
    return cond,val,maxdist

  end

end



























































#############################################################################

end
