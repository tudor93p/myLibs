module BandStructure
#############################################################################

using Distributed

import ..DlmF, ..LA, ..SpA, Arpack

import ..Algebra, ..Utils, ..ReadWrite



#===========================================================================#
#
# Decide between dense and sparse diagonalization
#
#---------------------------------------------------------------------------#

function get_eigen(H::Function, 
									 evect::Bool=false; 
									 tol::Float64=1e-8, 
									 nr_bands=nothing, 
									 sigma::Float64=tol/10)::Function 


  if isnothing(nr_bands)

    !evect && return k -> (LA.eigvals(LA.Hermitian(Array(H(k)))),nothing)
 
    
    return k-> LA.eigen(LA.Hermitian(Array(H(k)))) |> e->(e.values,e.vectors)
 
  end

	

  Areigs(k) = Arpack.eigs(H(k), nev=nr_bands, 
																sigma=Utils.Assign_Value(sigma, tol/10),
																tol=tol,
																ritzvec=evect)


	evect && return k -> Areigs(k) |> function (e)


													order = sortperm(real(e[1]))

#													println("\n",real(e[1])[order],"\n")
														
													return (e[1][order], e[2][:,order])


												end



	return k -> (sort(Areigs(k)[1],by=real),nothing)


end





#===========================================================================#
#
# Evals and operators on evects
#
#---------------------------------------------------------------------------#

function get_nrbands(nr_bands::Nothing, args...)::Nothing

	nothing 

end 


function get_nrbands(nr_bands::Int, H0::AbstractMatrix, lim::Float64=0.2)::Union{Nothing, Int}

	nr_bands < size(H0,1)*lim && return nr_bands

  return nothing

end



function Diagonalize(H, kPoints, filename=nothing; 
										 kLabels=nothing,
										 kTicks=[0],
										 filemethod="new", 
										 parallel=false, operators=[[],[]],
										 storemethod="dat",
										 dim=1,
										 tol=1e-8,  nr_bands=nothing, sigma=tol/10, kwargs...)

	k1 = selectdim(kPoints, dim, 1)

	rest_ks = Base.Iterators.drop(eachslice(kPoints, dims=dim), 1)

	nr_bands_calc = get_nrbands(nr_bands, H(k1))

  eig = get_eigen(H, !isempty(operators[1]); 
									tol=tol, nr_bands=nr_bands_calc, sigma=sigma)


	function restrict((e,p))

						# if no restriction desired 
						# if the restriction was implemented already 

		isnothing(nr_bands) | !isnothing(nr_bands_calc)	&& return (e,p)

						# if only the eigenvalues are calculated
						# if the   

		isnothing(p) | (nr_bands >= length(e)) && return (e, p)

		m,M = argmin(abs.(e .- sigma)) .+ [-1,1]*div(nr_bands,2)

		inds = max(m-1,0) + min(length(e)-M,0) .+ (1:nr_bands)
	
		return (e[inds], p[:,inds] ) # vectors come on columns 

	end

	out(k) = eig(k) |> restrict |> function	((e,p),)

															(e, [Op(p,k=k) for Op in operators[2]])

																end


	E1, Op1 = out(k1)
															
	bounds = cumsum([0;1;size.(Op1,2)])


  result = vcat(

		hcat(E1, Op1...),

		(parallel ? pmap : map)(rest_ks) do k

			out(k) |> out_k -> hcat(out_k[1], out_k[2]...)

    end...)




	Write!, outdict = ReadWrite.Write_NamesVals(filename,  storemethod, 
																		["Energy";operators[1]], result, bounds;
																		tol=tol, filemethod=filemethod)





	if size(kPoints,dim)==1
	
		return Write!("kLabels", reshape(range(0,1,length=size(result,1)),:,1), outdict)


	end	


	kLabels = Utils.Assign_Value(kLabels, range(0,1,length=size(kPoints,dim)))

	kLabels = reshape(repeat(kLabels,
													 inner=div(size(result,1), length(kLabels))),:,1)


	# kPoints should also be repeated

	Write!("kLabels", kLabels, outdict)

	Write!("kTicks", kTicks, outdict)


  return outdict


end



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#



function Read_Bands(filename,colors=[])

  return Utils.Read_NamesVals(filename,[["Energy","kLabels"];colors])

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

function occupied(E,P)

  return P[:,sortperm(E) .<= div(size(E,1),2)]

end

function WLO_(psi,kPoints)

  Psi = map(psi,eachrow(kPoints[1:end-1,:]))

  return Psi[1]'*mapreduce(P->P*P',*,Psi[2:end])*Psi[1]

end

function WLO(H,kPoints;tol=1e-8,nr_bands=nothing,sigma=tol/10)

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


