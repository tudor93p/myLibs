module BandStructure
#############################################################################

using Distributed

import ..DlmF, ..LA, ..SpA, Arpack

import ..Algebra, ..Utils, ..ReadWrite, ..Operators



#===========================================================================#
#
# Decide between dense and sparse diagonalization
#
#---------------------------------------------------------------------------#

function restrict_energies_inds(e::AbstractVector{Float64},
											 nr_bands::Int,
											 sigma::Float64)::Vector{Int}

	m,M = argmin(abs.(e .- sigma)) .+ [-1,1]*div(nr_bands,2)

	return max(m-1,0) + min(length(e)-M,0) .+ (1:nr_bands)

end  


function get_eigen(Hk::AbstractMatrix{<:Number},
									 sparsediag::Val{false},
									 evect::Val{false},
									 nr_bands::Nothing=nothing;
									 kwargs...
									 )::Tuple{Vector{Float64},Nothing}
	
	(LA.eigvals(LA.Hermitian(Matrix(Hk))), nothing)

end 

function get_eigen(Hk::AbstractMatrix{<:Number},
									 sparsediag::Val{false},
									 evect::Val{false},
									 nr_bands::Int;
									 sigma::Float64,
									 kwargs...
									 )::Tuple{Vector{Float64},Nothing}


	(e, p) = get_eigen(Hk, sparsediag, evect)

	return (e[restrict_energies_inds(e, nr_bands, sigma)], p)

end  



function get_eigen(Hk::AbstractMatrix{<:Number},
									 sparsediag::Val{false},
									 evect::Val{true},
									 nr_bands::Nothing=nothing;
									 kwargs...
									 )::Tuple{Vector{Float64},Matrix{ComplexF64}}
	
	eigen = LA.eigen(LA.Hermitian(Matrix(Hk)))

	return (eigen.values, eigen.vectors)

end 



function get_eigen(Hk::AbstractMatrix{<:Number},
									 sparsediag::Val{false},
									 evect::Val{true},
									 nr_bands::Int;
									 sigma::Float64,
									 )::Tuple{Vector{Float64},Matrix{ComplexF64}}

	e,p = get_eigen(Hk, sparsediag, evect) 

	nr_bands>=length(e) && return (e, p)

	inds = restrict_energies_inds(e, nr_bands, sigma) 

	return (e[inds], p[:,inds]) # vectors come on columns 

end 




function get_eigen(Hk::AbstractMatrix{<:Number},
									 sparsediag::Val{true},
									 evect::Val{false},
									 nr_bands::Int;
									 tol::Float64,
									 sigma::Float64,
									 )::Tuple{Vector{ComplexF64},Nothing}
  
	eigs = Arpack.eigs(Hk, nev=nr_bands, tol=tol, ritzvec=false,
										 sigma=Utils.Assign_Value(sigma, tol/10))

	return (sort(eigs[1], by=real), nothing)

end 


function get_eigen(Hk::AbstractMatrix{<:Number},
									 sparsediag::Val{true},
									 evect::Val{true},
									 nr_bands::Int;
									 tol::Float64,
									 sigma::Float64,
									 )::Tuple{Vector{ComplexF64}, Matrix{ComplexF64}}
  
	eigs = Arpack.eigs(Hk, nev=nr_bands, tol=tol, ritzvec=true, sigma=sigma)

	order = sortperm(eigs[1], by=real)

	return (eigs[1][order], eigs[2][:,order])

end 


function get_eigen(sparsediag::Val, args...; kw1...)::Function 

	function get_eigen_(Hk::AbstractMatrix; kw2...)::Tuple{Vector,Any}

		get_eigen(Hk, sparsediag, args...; kw1..., kw2...)

	end 

end 


function get_eigen(H::Function, args...; kwargs...)::Function 

	get_eigen(args...; kwargs...) ∘ H 

end 


#===========================================================================#
#
# Evals and operators on evects
#
#---------------------------------------------------------------------------#

function get_nrbands(H0::AbstractMatrix, lim::Float64=0.2;
										nr_bands=nothing, kwargs...)::Union{Nothing, Int}

	isa(nr_bands,Int) && nr_bands<size(H0,1)*lim && return nr_bands

  return nothing

end



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function prep_kLabels(kPoints::AbstractMatrix, dim::Int;
											kLabels=nothing, kwargs...)::AbstractVector

	if isnothing(kLabels)

		size(kPoints,dim)==1 && return [0]

		return range(0,1,length=size(kPoints,dim))

	elseif kLabels isa AbstractVector 

		@assert length(kLabels)==size(kPoints,dim) string(length(kLabels)," ",size(kPoints))

		return kLabels 

	end 

end 



function prep_operators(;operators=nothing, kwargs...
											 )::Tuple{Vector{String}, Vector} 

	(isnothing(operators) || isempty(operators)) && return (String[],Function[])

	@assert length(operators)==2 
	
	N,O = operators 
	
	@assert length(N)==length(O)

	return (N,O)


end  


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function Diagonalize_oneK(k::AbstractVector{Float64},
													Hk::AbstractMatrix{<:Number},
													lab,
													operators::Tuple{Vector,Vector},
													args...;
													sparsediag::Bool=false,
													kwargs...
													)::Dict{String,Any}


	e, p = get_eigen(Hk, Val(sparsediag), Val(!isempty(operators[1])), 
									 args...; kwargs...)

	return Dict{String,Any}("Energy"=>e,
													"kLabels"=>fill(lab, length(e)),
													(N=>Op(p; OpArg=k) for (N,Op)=zip(operators...))...)

end  


function Diagonalize_oneK(k::AbstractVector{Float64}, H::Function,
													args...; kwargs...)::Dict{String,Any}

	Diagonalize_oneK(k, H(k), args...; kwargs...)

end 


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function Diagonalize(H::Function, 
										 kPoints::AbstractMatrix{Float64}, 
										 filename::Nothing=nothing; 
										 filemethod::AbstractString="new", 
										 parallel::Bool=false, 
										 storemethod::AbstractString="jld",
										 dim::Int,
										 tol::Float64=1e-8,
										 sigma::Float64=tol/10, kwargs...)
	
	operators = prep_operators(;kwargs...)  

	kLabels = prep_kLabels(kPoints, dim; kwargs...)


	k1 = selectdim(kPoints, dim, 1)

	rest_ks = Base.Iterators.drop(eachslice(kPoints, dims=dim), 1)

	klab1, rest_klab = kLabels[1], kLabels[2:end]
	
	
	H1 = H(k1)

	nr_bands = get_nrbands(H1; kwargs...)

	DoK(k,h,l) = Diagonalize_oneK(k,h,l, operators, nr_bands; tol=tol, sigma=sigma)
	DoK((k,l)) = DoK(k,H,l)


	D = Utils.RecursiveMerge(
													 
				DoK(k1, H1, klab1),

				Utils.RecursiveMerge(DoK, zip(rest_ks, rest_klab); dim=dim);

				dim=dim, parallel=parallel)

	size(kPoints,dim)==1 && return D 

	for (N,V) in D 

		s = size(V)

		@assert length(s)>=2 "Strange dimension"

		S = [(s[1]*s[2],s[3:end]...), (s[1:end-2]..., s[end-1]*s[end])]

		D[N] = reshape(V, S[dim])#filter(>(1), S[dim]))

	end 

	haskey(kwargs, :kTicks) && setindex!(D, kwargs[:kTicks], "kTicks")

	return D

end


function Diagonalize(H::Function, 
										 kPoints::AbstractMatrix{Float64}, 
										 filename::Function;
										 storemethod::AbstractString="jld",
										 kwargs...)

	ReadWrite.Write_PhysObs(filename, storemethod, 
													Diagonalize(H, kPoints; kwargs...))

end 


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function prep_oper_names_rw(operators::Nothing=nothing, args...)::Vector{String}
	
	prep_oper_names_rw([], args...)

end 

function prep_oper_names_rw(operators::AbstractString,args...)::Vector{String} 
	prep_oper_names_rw([operators],args...) 

end 

function prep_oper_names_rw(operators::AbstractVector{<:AbstractString},
													 )::Vector{String} 

	operators 
end 

function prep_oper_names_rw(operators::AbstractVector{<:AbstractString},
														extra_opers
													 )::Vector{String}  

	[operators; prep_oper_names_rw(extra_opers)]

end 




function Read_Bands(filename::Function, operators, args...)::Dict 

	names = prep_oper_names_rw(operators, ["kLabels", "kTicks", "Energy"])
	
	return ReadWrite.Read_NamesVals(filename, names, args...)

end


function FoundFiles_Bands(filename::Function, operators, args...)::Bool

	names = prep_oper_names_rw(operators, ["kLabels","Energy"])

	return ReadWrite.FoundFiles_NamesVals(filename, names, args...)

end 
#===========================================================================#
#
# Time evolution of a state. psi(t) = e^iHt psi0
#
#---------------------------------------------------------------------------#

function TimeEvolution_expH(psi0,es,vs,ts,filename,operators=[[],[]];tol=1e-7)

  V = reshape(vs'*psi0,1,:).*vs

  states = V * exp.(1im*es.*reshape(ts,1,:))

  states[:,1] = psi0  

  weights(P;kwargs...) =  abs.(transpose(P)*conj(vs))



  Write, = ReadWrite.Write_NamesVals(filename; tol=tol)

  Write("Time",ts)
  Write("Eigenvals",es)

  for (name,oper) in zip([operators[1];"Coef"],[operators[2];weights])

    Write(name,oper(states))

  end


end




#===========================================================================#
#
# Analyzes the energy spectrum find a reasonable energy step
#
#---------------------------------------------------------------------------#


function Analyze_EnSplitting(En,kPoints)


  En = reshape(En,:,size(kPoints,1))

#  En[band_index,k_index]
# Algebra.OuterDist(E,E) |> d -> d[LA.triu(trues(size(d)),1)] |> d -> sort(d[d .> tol])


  bandsep = abs.(diff(En,dims=1))[:] #all band splittings

#  mesh = sort(vcat(meshsplit.(eachrow(En))...)) # all splittings En(k)->En(k')



  if size(kPoints,1) ==1

    closestk = 1

  else

    closestk = hcat(sortperm.(eachcol(Algebra.OuterDist(kPoints,kPoints)))...)[2,:]
  end

  

  p=2 

  D(i) = minimum(abs.(En[i:i,:] .- En[i-p.<=axes(En,1).<=i+p,closestk]),dims=1)
     #  splittings of level i|En[i,k] - En[closest i',closest k]| for all ks

  meshsplit = hcat(D.(axes(En,1))...)[:]# |> x -> sum(x)/length(x)

#  sort(x)[div(length(x),10)]

  return [sum(x)/length(x) for x in [meshsplit,bandsep]]
# DOS: delta should be smaller than the typical band separation (bandsep), but larger than the energy splitting due to the finite k-mesh (meshsplit)

end


#===========================================================================#
#
# Read the DOS
#
#---------------------------------------------------------------------------#


function Read_DOS(filename;weights="Lorentzian",LDOS=false)#::String)

  En,kPts = DlmF.readdlm.(filename.(["Energy","kPoints"]))

  scales = vcat(minimum(En),maximum(En),Analyze_EnSplitting(En,kPts)...)

  LDOS = LDOS && DlmF.readdlm(filename("LDOS"))

  return scales,fDOS(En,LDOS,weights)

end



function fDOS(En,LDOS=false;weights="Lorentzian",get_weights=false)

  weights = Algebra.get_Distribution(weights)


  function f1(E,delta)

    W = weights(E,En,delta)

    return sum(W,dims=2)

  end


  function f2(E,delta)

    W = weights(E,En,delta)
#    W = weights(E,transpose(En),delta)



    return W*LDOS./(sum(W,dims=2) .+ 1e-20)

  end

  function f3(E,delta)

    W = weights(E,En,delta)

    return sum(W,dims=2),W

  end


  function f4(E,delta)

    W = weights(E,En,delta)
#    W = weights(E,transpose(En),delta)

    return W*LDOS./(sum(W,dims=2) .+ 1e-20),W

  end


  if !get_weights

    !isa(LDOS,AbstractArray) && return f1

    return f1,f2

  else

    !isa(LDOS,AbstractArray) && return f3

    return f3,f4

  end

end













#===========================================================================#
#
# Wilson loop operators
#
#---------------------------------------------------------------------------#


#function WLO(H,kPoints,filename=nothing;kLabels=axes(kPoints,1),parallel=false,tol=1e-8,operators=[[],[]],nr_bands=nothing,sigma=tol/10)

function occupied(E::AbstractVector{Float64},
									P::AbstractMatrix{ComplexF64}; 
									dim::Int)::AbstractMatrix{ComplexF64}

	selectdim(P, [2,1][dim], sortperm(E) .<= div(length(E), 2) )

end

function WLO_(psi::Function, kPoints::AbstractMatrix{Float64}; dim::Int)::Matrix{ComplexF64}

	Psi = [psi(selectdim(kPoints, dim, i)) for i=1:size(kPoints,dim)-1]

  return Psi[1]'*mapreduce(P->P*P',*,Psi[2:end])*Psi[1]

end



function WLO(H::Function,
						 kPoints::AbstractMatrix{<:Float64};
						 tol::Float64=1e-8,
						 nr_bands=nothing,
						 sigma::Float64=tol/10)::Matrix{ComplexF64}

  eig = get_eigen(H,true;tol=tol,nr_bands=nr_bands,sigma=sigma)

  return WLO_(k->occupied(eig(k)...),kPoints)

end


function WLO_Spectrum(W,kPoints,filename=nothing;kLabels=axes(kPoints,1),parallel=false)

  evals(k) = angle.(LA.eigvals(W(k)))/2pi

  E = vcat((parallel ? pmap : map)(evals,eachrow(kPoints))...)


  function write_(name,matrix)

    
    if !isnothing(filename)
      open(filename(name),"w") do fout

        DlmF.writedlm(fout,ndims(matrix) == 1 ? reshape(matrix,:,1) : matrix)

      end
    end 

    return matrix
  end


  outdict = Dict()

  outdict["Energy"] = write_("Energy",E)

  labels = repeat(kLabels,inner=div(size(E,1),length(kLabels)))

  outdict["kLabels"] = write_("kLabels",labels)

  outdict["kPoints"] = write_("kPoints",kPoints)

  return outdict




end





function WLO2(H,W,kPoints;tol=1e-8,nr_bands=nothing,sigma=tol/10)

  eigH = get_eigen(H,true;tol=tol,nr_bands=nr_bands,sigma=sigma)

  eigWLO(k) = LA.eigen(W(k)) |> e-> (angle.(e.values),e.vectors)

  w(k) = occupied(eigH(k)...)*occupied(eigWLO(k)...)

  return WLO_(w,kPoints)

end




#===========================================================================#
#
# 
#
#---------------------------------------------------------------------------#























#############################################################################

end


