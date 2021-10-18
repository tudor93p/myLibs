module ReadWrite
#############################################################################

import ..DlmF
import FileIO, JLD

import ..Utils, ..ArrayOps


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


function Extension_Storemethod(storemethod::AbstractString)::String

	storemethod in ["dat","jld"] && return ".$storemethod"

	error("The 'storemethod' $storemethod is not supported")

end


function isLegendFile(storemethod::AbstractString)::Function

	storemethod == "dat" && return x->occursin("_Legend",x)

	storemethod == "jld" && return x->false

	error("The 'storemethod' $storemethod is not supported")

end 

function Name_fromLegend(storemethod::AbstractString)::Function

	storemethod == "dat" && return x->split(x,"_Legend")[1]
	
	storemethod == "jld" && return identity

	error("The 'storemethod' $storemethod is not supported")

end 

function LegendFile_fromName(storemethod::AbstractString)::Function

	storemethod == "dat" && return x->x*"_Legend"

	storemethod == "jld" && return identity

	error("The 'storemethod' $storemethod is not supported")

end 







#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#



function Write_NamesVals(filename::Function, 
												 storemethod::AbstractString="jld"; 
												 tol::Float64=1e-7, 
												 filemethod::AbstractString="new"
												 )::Tuple{Function,Dict}


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

	writable(s::AbstractString, name)::String = s


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
#				FileIO.save(FileIO.File(FileIO.format"JLD",fn), name, data) 

				FileIO.save(FileIO.File{FileIO.format"JLD"}(fn), name, data) 

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



function Read_NamesVals(filename::Function, Names::AbstractString, storemethod::AbstractString)::Dict

	Read_NamesVals(filename, [Names], storemethod)

end 


function good_data(D::T)::T where T<:AbstractDict 

	for v in values(D)
		
		v isa Union{<:Number, <:AbstractDict,<:AbstractArray} && continue 

		error("JLD data is not read as a dict! 'import JLD' in the script using ReadWrite should solve the problem")

	end 

	return D 

end 


function unique_existent_fileNames(filename::Function, Names::AbstractVector{<:AbstractString}, storemethod::AbstractString)::Tuple{Vector{String},Vector{String}}

	ufn,iufn = Utils.Unique(
							filename.(Names) .* Extension_Storemethod(storemethod);
													inds=:first)

	exist = isfile.(ufn) 
	
	return (Names[iufn[exist]], ufn[exist])

end 



function Read_NamesVals(filename::Function, Names::AbstractVector{<:AbstractString}, storemethod::AbstractString)::Dict{String,Any}
	
	n,fn = unique_existent_fileNames(filename, Names, storemethod)

	isempty(n) && return Dict{String,Any}()

	return good_data(Read_NamesVals_(n, fn, storemethod))

end 



function Read_NamesVals_(Names::T, fileNames::T, 
												 storemethod::AbstractString
												 )::Dict{String,Any} where T<:AbstractVector{<:AbstractString}

	if storemethod=="jld" 
		
		merge!(Dict{String,Any}(), JLD.load.(fileNames)...)

	elseif storemethod=="dat" 
		
		Dict{String,Any}(zip(Names, DlmF.readdlm.(fileNames)))

	end 

end




#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function FoundFiles_NamesVals(filename::Function, Names, 
															storemethod::AbstractString)::Bool
	
	fn = filename.(vcat(Names)).*Extension_Storemethod(storemethod)

	return all(isfile, unique(fn))

end




#===========================================================================#
#
# Prepare dict for writing by flattening its values 
# 										(if array, use cartesian inds for the new keys)
#
#---------------------------------------------------------------------------#

function flattenDict_AssertType(D::AbstractDict{K,V})::Int where {K,V}

	for (i,T) in enumerate([Number, AbstractArray])#, AbstractDict])

		(V<:T || isa(first(D).second,T)) && return i

	end 

	error()

	return 0

end 


function flattenDict_Iter(D::AbstractDict)

	T = flattenDict_AssertType(D) 

	T==1 && return (nothing,nothing,D,nothing)

	newkey(k,CIs...) = join([k;map(I->join(Tuple(I),"-"), CIs)...], "_")  

	T==2 && return ((k,I,V[I],newkey(k,I)) for (k,V) in pairs(D) 
																					for I in CartesianIndices(V))

	error()


end 


function flattenDict_Keys(D::AbstractDict)::Vector{<:Any}

	any(occursin.("_",keys(D))) && error("The keys already contain the character '_'. The result will not be read correctly.")

	T = flattenDict_AssertType(D) 

	T==1 && return collect(keys(D))

	T==2 && return [newk for (_,_,_,newk) in flattenDict_Iter(D)]

end 


function flattenDict(D::AbstractDict)::AbstractDict{<:Any, <:Number}

	T = flattenDict_AssertType(D)  

	T==1 && return D 
	
	T==2 && return Dict(k=>v for (_,_,v,k) in flattenDict_Iter(D))

end 

function flattenDict_Vals(D::AbstractDict, vals::Bool=true)::Vector{<:Number}

	vals || error()

	T = flattenDict_AssertType(D)  

	T==1 && return collect(values(D))

	T==2 && return [v for (_,_,v,_) in flattenDict_Iter(D)]

	error()

end 


function flattenDict_Vals(D::AbstractDict, K::AbstractString)::Number

	for (_,_,v,k) in flattenDict_Iter(D) 

		k==K && return v

	end

	error()

end 


function flattenDict_Vals(D::AbstractDict, Ks::Utils.List)::Vector{<:Number}

	flattenDict_AssertType(D)  

	newD = flattenDict(D) 

	return [newD[k] for k in Ks]
	
end								





#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function restoreDictEntries(matrix::AbstractMatrix, legend)::Dict{String, <:Any}
											
	K = string.(legend[:])

	!all(occursin.("_",K)) && return Dict(zip(K,matrix[:]))


	all_labels, all_inds = Utils.zipmap(legend) do L 

		left,right = split(L,"_")

		return (left, parse.(Int, split(right,"-")))
	
	end 


	return Dict(map(Utils.EnumUnique(all_labels)) do (label, loc)

		Pair(string(label),
				 ArrayOps.Array_from_ListIndsVals(all_inds[loc], matrix[loc]))
	end)


end 







#===========================================================================#
#
#	Write physical observable. Possibly with a legend file
#
#---------------------------------------------------------------------------#


function Write_PhysObs(filename::Function, storemethod::AbstractString;
											 kwargs...)

	get_legend(vals::AbstractDict)::Vector = sort(flattenDict_Keys(vals), by=string)
	get_legend(vals)::Nothing = nothing	 


	get_matrix(vals::AbstractVecOrMat{<:Number}, a...)::AbstractVecOrMat{<:Number} = vals

	get_matrix(vals::Number, a...)::Vector{<:Number} = [vals] 

	function get_matrix(vals::AbstractDict, legend::AbstractVector
											)::Vector{<:Number} 
		
		flattenDict_Vals(vals, legend)

	end 





	Write!, = Write_NamesVals(filename, storemethod; kwargs...)

	legend_file = LegendFile_fromName(storemethod)



	function new_item!((obs, vals), outdict::AbstractDict=Dict())

		new_item!(obs, vals, outdict)

	end



	function new_item!(obs, vals, outdict=Dict())
		
	
#		legend,matrix,data = legend_matrix_data(vals)
	
		if storemethod=="jld"

			Write!(obs, vals, outdict)

		elseif storemethod=="dat" 
			
			legend = get_legend(vals) 

			Write!(legend_file(obs), legend, outdict)


			Write!(obs, get_matrix(vals, legend)) 
								# ! careful ! this object is written, not returned 

			outdict[obs] = vals 
								# the original data 

		end

		return outdict

	end
	
	
	return new_item!,Dict()

end



function Write_PhysObs(filename::Function, storemethod::AbstractString,
											 target_dict::AbstractDict, 
											 target_obs=keys(target_dict);
											 kwargs...)

	Write!,outdict = Write_PhysObs(filename, storemethod)
	
	for obs in target_obs 
		
		Write!(obs, target_dict[obs], outdict)

	end

	return outdict

end 




#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

FoundFiles_PhysObs = FoundFiles_NamesVals


function Read_PhysObs(filename::Function, Names, storemethod::AbstractString)

	if storemethod=="jld"

		return Read_NamesVals(filename, Names, storemethod)

	elseif storemethod=="dat"

		legf = LegendFile_fromName(storemethod)

		obsf = Name_fromLegend(storemethod)

		out = Read_NamesVals(filename, [Names; legf.(vcat(Names))], storemethod)
													
		for legend in filter(isLegendFile(storemethod), keys(out))

			obs = obsf(legend)

			haskey(out, obs) || error("Legend exists, but no obs?")

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
