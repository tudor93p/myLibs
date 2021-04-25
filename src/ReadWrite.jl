module ReadWrite
#############################################################################

import ..DlmF
import FileIO#,JLD

import ..Utils


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function has_significant_imagpart(matrix; atol=1e-8, rtol=1e-5, mute=false)

  re_matrix = real(matrix)

	inds =  findall(abs.(re_matrix) .> atol)

	if !isempty(inds)

		rel_im_matrix = imag(Array(matrix)[inds])./Array(re_matrix)[inds]

		inds2 = findall(abs.(rel_im_matrix) .> rtol)
	
		if !isempty(inds2) 
		
			

			!mute && println("Number of elements with significant imaginary part: ",
																length(matrix[inds][inds2])
											)

			!mute && println("Maximum relative imaginary part: ",
																maximum(abs.(rel_im_matrix[inds2]))
											)

			return true 
		end
	end

	return false

end


#===========================================================================#
#
# Write results in different files
#
#---------------------------------------------------------------------------#

function Write_NamesVals(filename,  storemethod, Names, result, bounds; tol=1e-7, filemethod="new")


	Write!, outdict = Write_NamesVals(filename, storemethod;
																		filemethod=filemethod,
																		tol=tol)


  for (i,name) in enumerate(Names)

 		Write!(name, result[:,bounds[i]+1:bounds[i+1]], outdict)

  end

  return Write!, outdict

#  outdict = new_item("kLabels",reshape(repeat(kLabels,inner=div(size(result,1),length(kLabels))),:,1),outdict)

#  outdict = new_item("kPoints",kPoints,outdict)

#  return outdict


end





#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function Extension_Storemethod(storemethod)

	storemethod in ["dat","jld"] && return ".$storemethod"

	error("The 'storemethod' $storemethod is not supported")

end


function isLegendFile(storemethod)

	storemethod == "dat" && return x->occursin("_Legend",x)

	storemethod == "jld" && return x->false

	error("The 'storemethod' $storemethod is not supported")

end 

function Name_fromLegend(storemethod)

	storemethod == "dat" && return x->split(x,"_Legend")[1]
	
	storemethod == "jld" && return identity

	error("The 'storemethod' $storemethod is not supported")

end 

function LegendFile_fromName(storemethod)

	storemethod == "dat" && return x->x*"_Legend"

	storemethod == "jld" && return identity

	error("The 'storemethod' $storemethod is not supported")

end 







#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#



function Write_NamesVals(filename, storemethod="jld"; 
												 tol=1e-7, filemethod="new")


	function writable(matrix::AbstractArray{<:Number}, name)

		if has_significant_imagpart(matrix; atol=tol)
			
			println("Observable '$name'")

			!isnothing(filename) && println(filename(name))

			println()

			@warn "Observable '$name' not real!"

		end

		return real(matrix)

	end

	writable(matrix::AbstractArray{<:AbstractString},name) = matrix

	writable(matrix::AbstractArray{<:Char},name) = matrix

#	writable(matrix::Char, name) = matrix

	writable(x::Number, name) = writable(vcat(x), name) 

	function writable(D::AbstractDict, name) 
		
		Dict(k=>writable(v,name) for (k,v) in pairs(D))

	end 



	ext = Extension_Storemethod(storemethod)


  function write_(name, data)

    if !isnothing(filename) 

			fn = filename(name)*ext

			if storemethod=="dat"

	      open(fn, filemethod=="append" ? "a" : "w") do fout

  	       DlmF.writedlm(fout, data)

    	  end

			elseif storemethod=="jld"

#				JLD.save(FileIO.File(FileIO.format"JLD",fn), name, data)
				FileIO.save(FileIO.File(FileIO.format"JLD",fn), name, data)

			end 

    end
    
    return data 
		
  end



	function new_item!(name::AbstractString, val, outdict::AbstractDict=Dict())

#		@show name typeof(val) keys(outdict)

		if !isnothing(name) & !isnothing(val)

			outdict[name] = write_(name, writable(val,name))

		end 

    return outdict

  end

	function new_item!((name, val)::T, outdict::AbstractDict=Dict()) where T<:Tuple{AbstractString,Any}

		new_item!(name, val, outdict)

	end 	

  return new_item!,Dict()

end





#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function Delete_NamesVals(filename, Names, storemethod)

	ext = Extension_Storemethod(storemethod)

	get_legend = LegendFile_fromName(storemethod)

	Ns = isa(Names,AbstractString) ? [Names] : Names

#	println("\nDelete ",join(Ns,", "),"\n")

	for fn in filename.(Ns)

		rm.(filter(isfile, unique([fn,get_legend(fn)]).*ext))

	end 

end






#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#



function Read_NamesVals(filename, Names, storemethod)
	
	Names = isa(Names,AbstractString) ? [Names] : Names

	fileNames = unique(filename.(Names).*Extension_Storemethod(storemethod))


	if storemethod=="jld"

		FNs = filter(isfile, fileNames)

		isempty(FNs) && return Dict()

		return merge(map(FNs) do fn 

				FileIO.load(FileIO.File(FileIO.format"JLD",fn))

		end...) 

	elseif storemethod=="dat"

	  outdict = Dict()
	
		for (n,fn) in zip(Names,fileNames)

			if isfile(fn)

				outdict[n] =DlmF.readdlm(fn)

			end 
	
	  end
	
	  return outdict

	end 

end




#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function FoundFiles_NamesVals(filename, Names, storemethod)
	
	fn = filename.(vcat(Names)).*Extension_Storemethod(storemethod)

	return all(isfile, unique(fn))

end




#===========================================================================#
#
# Prepare dict for writing by flattening its values 
# 										(if array, use cartesian inds for the new keys)
#
#---------------------------------------------------------------------------#

function flattenDictEntries(D::AbstractDict; onlykeys=false, 
																							onlyvals=nothing,
																							)

	any(occursin.("_",keys(D))) && error("The keys already contain the character '_'. The result will not be read correctly.")

	T = typeof( first(D).second )


	onlykeys && !isnothing(onlyvals) && error("The kwargs 'onlykeys' and 'onlyvals' are mutually exclusive")

#	getkeys(D) = [D[k] for k in Keys]

	if T<:Number 
		
		onlykeys && return collect(keys(D))
		
		isnothing(onlyvals) && return D

		onlyvals isa Bool && return collect(values(D))

		onlyvals isa AbstractString && return D[onlyvals]

		isList(onlyvals, AbstractString) && return [D[k] for k in onlyvals]

		error("'onlyvals=$onlyvals' not supported")

	end 
	
	T<:AbstractArray || error("$T not supported")

	key(k,I) = string(k, "_", join(Tuple(I),"-")) 

	iter = ((k,I,V[I]) for (k,V) in pairs(D) for I in CartesianIndices(V))


	onlykeys && return [key(k,I) for (k,I,v) in iter]

	onlyvals isa Bool && return [v for (k,I,v) in iter]

	if onlyvals isa AbstractString 
	
		for (k,I,v) in iter 

			key(k,I)==onlyvals && return v

		end

	end



	out = Dict(key(k,I)=>v for (k,I,v) in iter)

	isnothing(onlyvals) && return out

	isList(onlyvals, AbstractString) && return [out[k] for k in onlyvals]


	error("'onlyvals=$onlyvals' not supported")


end								



#function restoreDictEntries(d::AbstractDict)
#
#	!occursin("-",first(keys(d))) && return d 
#
#	K = collect(keys(d))
#
#	sK = hcat(split.(K,"-")...)
#
#	return Dict(map(unique(sK[1,:])) do label 
#
#		select = label.==sK[1,:]
#		
#		return label => Array_from_ListIndsVals(
#												transpose(parse.(Int, sK[2:end, select])),
#												[d[k] for k in K[select]])
#	end)
#
#end 

function restoreDictEntries(matrix::AbstractMatrix, legend)
											
	K = string.(legend[:])

	!all(occursin.("_",K)) && return Dict(zip(K,eachcol(matrix)))

	sK = hcat(map(legend) do L
					
						left,right = split(L,"_")

						return vcat(left,split(right,"-"))

					 end...)

#	column: [label, inds...]


	return Dict(map(unique(sK[1,:])) do label 

		select = label.==sK[1,:] 	# all elements corresp to label
		
		return label => Array_from_ListIndsVals(
												transpose(parse.(Int, sK[2:end, select])),
												eachcol(matrix[:,select])
																						)
	end)

end 







#===========================================================================#
#
#	Write physical observable. Possibly with a legend file
#
#---------------------------------------------------------------------------#


function Write_PhysObs(filename, storemethod; tol=1e-7)

	function assert_type(v)

		isa(v, Number) && return 1 
		
		isa(v, AbstractDict) && return 3 

		isList(v, Number) && return 2 


		error(typeof(v)," not supported")

	end 


	concat(f,vals) = vcat([reshape(f(v),1,:) for v in vals]...)
# above: preallocate and avoid vcat ?

	function get_legend(vals)

		assert_type(vals[1]) in [1,2] && return nothing
		
		return sort(flattenDictEntries(vals[1], onlykeys=true))

	end 

	function get_matrix(vals)

		assert_type(vals[1]) in [1,2] && return concat(vcat, vals)

		return concat(v->flattenDictEntries(v; onlyvals=get_legend(vals)), vals)


	end 

	function get_data(vals)

		assert_type(vals[1]) in [1,2] && return concat(vcat, vals)

		K = collect(keys(vals[1]))

		(T,S,I0) = vals[1][K[1]] |> a -> (eltype(a), size(a), fill(:, ndims(a)))

		data = Dict(k=>zeros(T, (length(vals), S...)) for k in K) 

		for (i,val) in enumerate(vals), (k,v) in pairs(val) 

			data[k][i, I0...] = v

		end

		return data

	end



	Write!, = Write_NamesVals(filename, storemethod; tol=tol)

	legend_file = LegendFile_fromName(storemethod)



	function new_item!((obs, vals), outdict::AbstractDict=Dict())

		new_item!(obs, vals, outdict)

	end



	function new_item!(obs, vals, outdict=Dict())
		
	
#		legend,matrix,data = legend_matrix_data(vals)
	
		if storemethod=="jld"

			#return 
			Write!(obs, get_data(vals), outdict)

		elseif storemethod=="dat"

			Write!(obs, get_matrix(vals)) 
								# ! careful ! this object is written, not returned 

			Write!(legend_file(obs), get_legend(vals), outdict)

			outdict[obs] = get_data(vals)

		end

		return outdict

	end
	
	
	return new_item!,Dict()

end




#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

FoundFiles_PhysObs = FoundFiles_NamesVals


function Read_PhysObs(filename, Names, storemethod)

	if storemethod=="jld"

		return Read_NamesVals(filename, Names, storemethod)

	elseif storemethod=="dat"

		legf = LegendFile_fromName(storemethod)

		obsf = Name_fromLegend(storemethod)

		out = Read_NamesVals(filename, [Names; legf.(vcat(Names))], storemethod)
													
		for legend in filter(isLegendFile(storemethod), keys(out))

			obs = obsf(legend)

			!haskey(out, obs) && error("Legend exists, but no obs?")

			out[obs] = restoreDictEntries(pop!(out, obs), pop!(out, legend))
	
		end
	
		return out 

	end 

end 





#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#





function ChangeStoreMethod_PhysObs(filename, Names, source, dest;
																	delete=true)


	if source!=dest 

		Write!, = Utils.Write_NamesVals(filename,  dest)

#		@show Names vcat(Names) source isLegendFile(source).(vcat(Names))


		for name in filter(!isLegendFile(source), vcat(Names))

#			println("\nMoving '$name' from '$source' to '$dest'\n")

#			@show 1,name 
#			M = Read_PhysObs(filename, name, source)[name]
#			@show 2,name
#			M isa AbstractDict && @show keys(M) size.(values(M))
#			M isa AbstractArray && @show size(M)

#			@show 3,name 
#			@show typeof(M)

#			Write!(name, M)
			Write!(name, Read_PhysObs(filename, name, source)[name])

			delete && Delete_NamesVals(filename, name, source)

		end 

	end

	return FoundFiles_PhysObs(filename, Names, dest)

#	return true # if succesful, files exist 

end 















































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































#############################################################################
end
