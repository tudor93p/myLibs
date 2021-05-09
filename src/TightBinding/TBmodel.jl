module TBmodel
#############################################################################


import ..LA, ..SpA

import ..Utils, ..Algebra


#===========================================================================#
#
# Calculates 	--> the Hamiltonian of the layer parallel to a given surface
#							--> the perpendicular hopping Hamiltonian
#
#---------------------------------------------------------------------------#

function BlochHamilt_Parallel(pyLatt,arg1_BH,arg2_BH,surface=0)
 
  Parallel_Layer = pyLatt.Reduce_Dim(surface)
      
      # -1 if "surface" is given for julia instead of python
      # 	(in julia it counts from 1, in python from 0)
      # Remove the specified dimension "surface" => dimension = d-1
  
  return Bloch_Hamilt(arg1_BH(Parallel_Layer);arg2_BH...,argH="phi")

  	  # a function of k parallel, vector of dimension d-1
  
end

 
function BlochHamilt_Perpendicular(pyLatt,arg1_BH,arg2_BH,surface=0)

  Coupled_Layers = pyLatt.Reduce_Dim(surface,complement=true)
        # Remove from Latt all dimensions not in "surface" 
				# the resulting lattice has dimension 1
				
  (ms,Rms,AtomsUC) = arg1_BH(Coupled_Layers)

	ms,Rms = map([ms,Rms]) do A

		return ndims(A) == 2 ? Utils.ArrayToList(A) : A
	end


	for sgn in [+1,-1]

		useful_ms = findall([m[1]*sgn > 0 for m in ms])

		isempty(useful_ms) && continue 
		
		return Dict(map(useful_ms) do i
	
			m,Rm = sgn*ms[i][1],hcat(sgn*Rms[i]...)

			return m => HoppingMatrix(AtomsUC, AtomsUC .+ Rm; arg2_BH...)

		end)

	end

	return Dict()

end




#===========================================================================#
#
#		Defines the arrangement of the Hamiltonian:
#
#					[ all orbitals of atom 1, all orbitals of atom 2, ... ]
#				
#---------------------------------------------------------------------------#

function Hamilt_indices(orbital,atom,nr_orbitals)

  return orbital .+ (atom .- 1) * nr_orbitals

end

function Hamilt_indices_all(orbs,atoms,d0=length(orbs);
														iter="atoms",flat=false)

  iter=="orbitals" && flat==true && println("\nWarning! The function returns flat  Hamiltonian indices by iterating over orbitals. In multi-orbital system, the order will *not* be consistent with the way the Hamiltonian is built and unexpected results will occur.\n")


  any(isnothing.([orbs,atoms,d0])) && error("Not enough info provided")

  f(x) = flat ? vcat(x...) : x

  iter == "orbitals"  && return f(map(i->Hamilt_indices(i,atoms,d0),orbs))
  
  iter == "atoms" && return f(map(i->Hamilt_indices(orbs,i,d0),atoms))

  error("'iter' should be either 'orbitals' or 'atoms'")

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


function flat_indsvals(i::Int64, j::Int64, v::AbstractMatrix, d0::Int64=1)::Matrix

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

function get_hopps(atom_i::Int, Ri::AbstractVector, Rsj::AbstractMatrix, hopping::Function, hopp_cutoff::Float64, d0::Int)::Matrix{ComplexF64}


  hopps = hopping.([Ri],eachcol(Rsj))
		# compute the hopping amplitudes from Ri to all Rj-s


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
											 Hopping::Function, Nr_Orbitals=1, 
											 parallel=false, hopp_cutoff=1e-6, kwargs...
											 )::SpA.SparseMatrixCSC

# Nr_Orbitals = trunc(Int64,sqrt(length(Hopping(Rsi[1,:],Rsj[1,:]))))

  indsvals_to_Tm(vcat(

      (parallel ? pmap : map)(enumerate(eachcol(Rsi))) do (i,Ri)
      	# for each atom i at position Ri

				return get_hopps(i, Ri, Rsj, Hopping, hopp_cutoff, Nr_Orbitals)

	  	# find the atoms where it hopps and the amplitudes

    end...),size(Rsi,2)*Nr_Orbitals, size(Rsj,2)*Nr_Orbitals)

end

#===========================================================================#
#
# 	Computes all tight-binding matrices
#				for a certain lattice and hopping function
#
#---------------------------------------------------------------------------#

function Compute_Hopping_Matrices(
						Lattice::NTuple{3,AbstractArray};
						dist_tol=1e-5,parallel=false,	kwargs...
						)::Tuple{SpA.SparseMatrixCSC, Any}


#  ms, Rms, AtomsUC = Lattice 


  # 'AtomsUC' 	-> the positions of the atoms in the central unit cell
  # 'Rms'	-> the positions of the unit cells around the central one
  # 'ms'	-> the integer vectors 𝐦  designating those unit cells
  	#	-> one 𝐦 vector on each column of the matrix ms
	#	-> 𝐦 has the same length as the number of lattice vectors
	# 	-> each UC in the list is 𝐑_𝐦 = Σ_j (m_j ⋅ 𝐚_j)

  # kwargs include
  # Hopping 	-> the function f(ri,rj) which gives the hopping amplitude
  #	 		where the atoms in the central cell can hop;
  # hopp_cutoff -> the cutoff for the hopping amplitudes
  # Nr_Orbitals -> the size of the matrix f(ri,rj) 
  # ...

	AtomsUC = Lattice[3] 

	ms,Rms = map(Lattice[1:2]) do A 

		ndims(A)<2 && return A 

		ndims(A)==2 && return collect(eachcol(A))

		error()

	end 

#  ms,Rms = map(A -> ndims(A) == 2 ? Utils.ArrayToList(A) : A, [ms,Rms])

	# if few procs, Rms will be distributed. Otherwise -> atoms
	
  Tms = (parallel && nprocs()<=size(Rms,2) ? pmap : map)(Rms) do Rm

    HoppingMatrix(AtomsUC,	Rm .+ AtomsUC;
													kwargs...,
      										parallel = parallel && nprocs()>size(Rms,2)
      										)
  end


  nonzeroT = findall(SpA.nnz.(Tms) .> 0)
		# check which T_𝐦  contain only zeros


	intra = findfirst(m -> m==-m, ms)
										
		# corresponds to the intra-cell hopping

  inter = nonzeroT[nonzeroT .!= intra]
		# correspond to the non-zero inter-cell hopping

  inter_out = !isempty(inter) ? map(A->A[inter], [ms,Rms,Tms]) : nothing

  return Tms[intra], inter_out

end


#===========================================================================#
#
# Compute the TB matrices, (save them) and return the Hamiltonian
#
#---------------------------------------------------------------------------#

function Bloch_Hamilt(Lattice; argH::String="",savefile="",	Hopping...)
			
  data = Compute_Hopping_Matrices(Lattice; Hopping...)

  isempty(savefile) || error("Save_Hamilt_Data not yet rewritten")

	# Save_Hamilt_Data(savefile,mRT...)

  return Assemble_Bloch_Hamiltonian(data...; argH=argH)

end




function Bloch_Hamilt(savefile::String="";argH::String="")

  error("Read_Hamilt_Data not yet rewritten")

  return Assemble_Bloch_Hamiltonian(Read_Hamilt_Data(savefile)...,argH=argH)

end


#===========================================================================#
#
# Build a tight-binding Hamiltonian H(k) based on hopping matrices
#
#---------------------------------------------------------------------------#


function Assemble_Bloch_Hamiltonian(intra::AbstractMatrix, inter::Nothing=nothing; kwargs...)::Function

	# if this is a finite system or the Hamiltonian is simply local
	
	H_no_k(args...) = intra 

end 
	
function Assemble_Bloch_Hamiltonian(intra::AbstractMatrix, inter::Utils.List;
																		argH::String="")::Function

  (ms,Rms,Tms) = inter

  argH = get_argH(argH, ms, Rms)

  iter = zip(Tms, [Rms,ms][findfirst(argH .== ["k","phi"])])


  return function h(k)
	
    h1 = sum([T*exp(Algebra.dot(k,R)im) for (T,R) in iter])

					# Bloch phase ℯ^(i 𝐤 ⋅ 𝐑 ), even if dimensions mismatch
   
    return intra + h1 + h1'

  end




end



#===========================================================================#
#
#		Helper function for deciding what argument the Hamiltonian gets
#
#---------------------------------------------------------------------------#

function get_argH(argH,ms,Rms)
  
  args = ["k","phi"]

  argH in args && return argH

	isempty(argH) || error("** Error in Bloch_Hamiltonian: argH=$argH not understood. It should be an empty string or either of ",join((a->"'$a'").(args),", ")," **")

	# it should be physical 'k' if the length of the lattice vectors,
	#	 is equal to their number, and Bloch phase 'phi' otherwise

  length(ms[1]) == length(Rms[1]) && return args[1]

  return args[2]

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



function Update_Maxdist(dist_tol,hopp_cutoff)


  function update(maxdist,param,d)

    abs(param) > hopp_cutoff && return max(maxdist,d)+2dist_tol
  
    return maxdist

  end

  return 2dist_tol,update

end


#===========================================================================#
#
# Initialize hopping terms -- different methods
# 		 
#
#---------------------------------------------------------------------------#


function Add_Hopping_Terms(d0,hopp_cutoff::Float64=1e-8)

  small = Utils.fSame(hopp_cutoff)

  matrixval(val) = Matrix{Complex{Float64}}(hcat(val))

  Zero = zeros(Complex{Float64},d0,d0)


  function Append(hoppings::Vector{Function},value,args...)::Vector{Function}

    small(value) && return hoppings
  

    return [hoppings; Hopping_Term(matrixval(value),Zero,args...)]

  end



  function Sum(hoppings::Vector{Function},condition)

    length(hoppings) == 0 && return (ri,rj) -> Zero

#    f(ri,rj) = reduce(+, map(h -> h(ri,rj), hoppings))

		f(ri,rj) = mapreduce(h->h(ri,rj), +, hoppings)

    return Hopping_Term(matrixval(1), Zero, condition, f)

  end



  return Append, Sum, Function[]


end


	# ------------ General case, matrix ------------ #

function Hopping_Term(value::AbstractMatrix,Zero::AbstractMatrix,condition::Function,fun::Function)

    return function f(ri::AbstractVector,rj::AbstractVector)

      condition(ri,rj) == true ? value .* fun(ri,rj) : Zero

    end

end

	# ------------ Constant matrix if condition  ------------ #

function Hopping_Term(value::AbstractMatrix,Zero::AbstractMatrix,condition::Function,fun::AbstractMatrix)

    return function f(ri::AbstractVector,rj::AbstractVector)

      condition(ri,rj) == true ? value .* fun : Zero

    end

end


	# ------------ Always return the hopping, function ------------ #

function Hopping_Term(value::AbstractMatrix,Zero::AbstractMatrix,condition::Bool,fun::Function)

  return function f(ri::AbstractVector,rj::AbstractVector)
  
    value .* fun(ri,rj)

  end

end

	# ------------ Always return the hopping, const ------------ #

function Hopping_Term(value::AbstractMatrix,Zero::AbstractMatrix,condition::Bool,fun::AbstractMatrix)

  return function f(ri::AbstractVector,rj::AbstractVector)
  
    value .* fun

  end

end

	# ------------ Constant hopping, matrix  ------------- #


function Hopping_Term(value::AbstractMatrix,Zero::AbstractMatrix,condition::Function)


  return function f(ri::AbstractVector,rj::AbstractVector)

    condition(ri,rj) == true ? value : Zero

  end

end



	# ------------ Hopping term from TB file ------------- #

function Hopping_Term_fromTBfile(filename)

  #TBmodel.Save_Hamilt_Data(filename,ms,Rms,Tms) #was done before

  Hopping_Term_fromTBmatr(Read_Hamilt_Data(filename)...)

end


function Hopping_Term_fromTBmatr(ms,Rs,Ts)

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
