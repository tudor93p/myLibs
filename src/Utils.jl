module Utils
#############################################################################



import ..LA 


import Dates, Combinatorics



using OrderedCollections:OrderedDict
import Random


const List = Union{AbstractVector, AbstractSet, Tuple, Base.Generator}




#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function Backtracking(data,
#											root::Function,
											possible_extensions::Function, 
											promising_candidate::Function, 
											accept_sol::Function,
											solutions::Vector{Dict}=Dict[],
											)::Vector{Dict}

	Backtracking(data,
#							 root,
							 possible_extensions, 
							 promising_candidate, 
							 accept_sol, 
							 (data, solutions, candidate) -> push!(solutions, candidate),
							 solutions
							)
end 


function Backtracking(data,
#											root::Function,
											possible_extensions::Function, 
											promising_candidate::Function, 
											accept_sol::Function, 
											output::Function,
											solutions::Vector{Dict}=Dict[],
											)::Vector{Dict}

	function backtrack!(data, solutions, candidate)


		promising_candidate(data, candidate) || return true 

		if accept_sol(data, candidate) 
			
			carry_on = output(data, solutions, candidate)
			
			return !isa(carry_on,Bool) || carry_on

		end 


		for extension in possible_extensions(data, candidate)

			carry_on = backtrack!(data, solutions, extension)

			isa(carry_on,Bool) && !carry_on && return false 
	
		end 

	end 


	backtrack!(data, solutions, Dict())#root(data))

	return solutions 

end 








#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function DistributeBallsToBoxes(balls::Int, boxes::Int)

	balls<0 && return -DistributeBallsToBoxes(-balls, boxes)

	map(Combinatorics.combinations(1:(balls+boxes-1), boxes-1)) do d

		diff(vcat(0, d, balls+boxes)) .- 1

	end

end 

function DistributeBallsToBoxes(list::Any, boxes::Int)

	DistributeBallsToBoxes(size(list,1), boxes)

end 

function EqualDistributeBallsToBoxes(balls::Int, boxes::Int)::Vector{Int}

  div(balls,boxes) .+ (1:boxes.<=balls%boxes)

end 

function EqualDistributeBallsToBoxes(list::Any, boxes::Int)::Vector{Int}

	EqualDistributeBallsToBoxes(size(list,1), boxes)

end 



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function ReplaceByNeighbor!(bad_element::Function, A::AbstractArray, start=CartesianIndex(first.(axes(A))))
	
	cart_ind = findnext(bad_element, A, CartesianIndex(start)) 

	isnothing(cart_ind) && return A


	for k in 1:maximum([sum(abs, cart_ind.I .- C) for C in [1, size(A)]])
					# neighbors 

		for J in sort(DistributeBallsToBoxes(k, ndims(A)), by=LA.norm)

			for S in Base.product(fill([1,-1],ndims(A))...)	# Signs

				new_cart_ind = CartesianIndex(cart_ind.I .+ Tuple(S.*J))

				all(1 .<= new_cart_ind.I .<= size(A)) || continue 

				a = A[new_cart_ind]

				bad_element(a) && continue

				A[cart_ind] = a 

				return ReplaceByNeighbor!(bad_element, A, cart_ind)

			end 

		end 

	end 
				
	error("The whole array has 'unwanted' elements")

end 


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#
is_float(S) = any(Ti -> S<:Ti, [Float64, Complex{<:Float64}])

is_exact(S) = any(Ti -> S<:Ti, [Integer, Rational, AbstractString, AbstractChar, Complex{<:Integer}, Complex{<:Rational}])



function Unique!(V::AbstractArray{T}; sorted=false, kwargs...) where T

	i = Unique(V; kwargs..., sorted=false, inds=:first)[2]

	deleteat!(V, filter(!in(i), axes(V,1)))

	sorted && sort!(V)

end

function tolNF(tol::Real)

	if isa(tol,Int)
		
		return tol,10.0^(-tol)

	elseif isa(tol,Float64)

		tol==0.0 && return 50, tol
	
		return Int(ceil(-log10(tol))), tol

	end 

end 


function Unique(V::AbstractArray{T}; 
								tol=1e-8, inds=nothing,
								sorted=false,
								check_type=true) where T



#	v isa AbstractArray &&
#
#	isList(v) || error("Type ",typeof(v)," not supported.")

	if !check_type || is_exact(T) || all(is_exact ∘ typeof, V)

		U = (sorted ? sort : identity)(unique(V))
	
		isnothing(inds) && return U

		i_f = findfirst(inds.==[:first,:all,:last])

		isnothing(i_f) && error("Cannot get indices with '$inds'")

		f = [findfirst, findall, findlast][i_f]

#		typeof(inds)<:Function || error("Cannot get indices with ", typeof(inds))
	
		return U, [f(isequal(u), V) for u in U]

	end 


	 
	function get_newV()

		ntol,ftol = tolNF(tol)

		ftol==0.0 && return V 

		is_float(T) && return trunc.(V, digits=ntol)

		fs, newT = typeof.(V) |> t -> (findall(is_float, t), promote_type(t...))

		W = AbstractArray{newT}(V)
	
		W[fs] = trunc.(W[fs], digits=ntol)

		return W

	end

	I = Unique(get_newV(); tol=tol, inds=Assign_Value(inds, :first),
						 						 sorted=sorted, check_type=false)[2]

	isnothing(inds) && return V[I]

	return V[first.(I)], I


end 

#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#



function flatmap(f,itr)

	vcat(map(f,itr)...)

end

function flatmapif(f, pred, itr)

	vcat(mapif(f, pred, itr)...)

end 

#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function mapif(f::Function,pred::Bool,itr)

	return pred ? map(f,itr) : []

end

function mapif(f::Function,pred::Function,itr)

	filter(pred,map(f,itr))

end

function zipmap(args...)

	zip(map(args...)...)

end 

function Zipmap(args...)

	Zip(map(args...)...)

end 

function filterzip(f, itr)

	filter(f, zip(itr...))

end 

function Zip(Args...) # zip for args of different lengths 

	map(1:maximum(length, Args)) do i

		[A[min(i,end)] for A in Args]

	end 


end 


#function zipifmap(args...)
#
#	zip(mapif(args...)...)
#
#end 
#
#
#function ifmapzip(f, pred, tr)
#
#	filter(pred, map(f, zip(itr)))
#
#end 
#
#function ifzipmap(f::Function, pred::Function, itr)
#	
#	filter(pred, zip(map(f,itr)...))
#
#end



invmap(arg, fs...) = invmap(arg, fs)

function invmap(args, fs)

	map(f->f(args...), fs)

end 









#===========================================================================#
#
#	Apply function on a given axis. Iterate the rest 
#
#---------------------------------------------------------------------------#


# Base function mapslices(f, y; dims=dim)


#===========================================================================#
#
#	Identify adjacent identical integers in a list
#
#---------------------------------------------------------------------------#

#[0,0,0,0,0,0,3,1,1,1] => [1:6,7:7,8:10]

function IdentifySectors(list,sectors=[],start=0; tol=1e-8)

	for i in axes(list,1)

		if isapprox(list[i],list[1],atol=tol)		# within the sector
		
			i==length(list) && return vcat(sectors,[1+start:i+start])
																			# sector ended abruptly

		else  # sector ended at i-1

			return IdentifySectors(list[i:end],
														 vcat(sectors,[1+start:i-1+start]),
														 start+i-1;
														 tol=tol)

		end
	end

end




#===========================================================================#
#
#	Identify Ranges -- group together inds which form ranges 
#
#---------------------------------------------------------------------------#


function IdentifyRanges(list)

	D2 = (diff(diff(list)).==0)

	!any(D2) && return list

	sectors =	map(filter(s->all(D2[s]),IdentifySectors(D2))) do s

							return minimum(s):maximum(s)+2
						end


	conflict(S) = (j -> maximum(S[j-1])>=minimum(S[j]))

	while true 

		i = findfirst(conflict(sectors),2:length(sectors))

		isnothing(i) && break

		s1,s2 = sectors[i:i+1]

		if length(s1) < length(s2)
			
			sectors[i] = minimum(s1):minimum(s2)-1

		else

			sectors[i+1] = maximum(s1)+1:maximum(s2)

		end

		sectors = filter(s->length(s)>2,sectors)

	end

#	if isempty(sectors)
#
#		d = filter(di->di>0,Algebra.FlatOuterDiff(q,q)[:])
#		
#		u = unique(d)
#		
#		c = zeros(Int64,length(u))
#		
#		for n in d c[indexin(n,u)[1]] +=1  end
#
#	end

	sector(i) = findfirst(s->i in s,sectors) 

	vcat(map(IdentifySectors(map(sector,axes(list,1)))) do S

		!in(S,sectors) && return list[S]

		(a,b) = extrema(S)
		
		step = list[a+1]-list[a]

		return [step==1 ? (list[a]:list[b]) : (list[a]:step:list[b])]

	end...)


end




#===========================================================================#
#
# check if arg is a vector or a tuple 
#
#---------------------------------------------------------------------------#




function isTuple(arg::T, Ti=Any) where T

	if T<:Type 

		arg <: Tuple{Vararg{<:Ti}} && return true 

	elseif T<:Tuple{Vararg{<:Ti}}
		
		return true 

	elseif T<:Tuple 

		for a in arg 

			typeof(a)<:Ti || return false 
			
		end 

		return true 

	end 

	return false 

end



function isList(arg::T, Ti=Any) where T

	for S in [AbstractVector, AbstractSet]

		if T<:Type 
			
			arg<:S{<:Ti} && return true 

		elseif T<:S{<:Ti}

			return true 

		elseif T<:S 

			for a in arg 

				typeof(a)<:Ti || return false 
				
			end 

			return true 

		end 

	end 


	return isTuple((typeof(arg) <: Base.Generator) ? Tuple(arg) : arg, Ti)

end


function isList(;T=Any)

	arg -> isList(arg, T)
	
end 







#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function flat(list_of_lists...; keep=nothing)

	i = findfirst(isList, list_of_lists)

	if isnothing(i) 
		
		isnothing(keep) && return vcat(list_of_lists...)
		
		return filter(keep, vcat(list_of_lists...))

	end

	return flat(list_of_lists[1:i-1]..., 
							list_of_lists[i]..., 
							list_of_lists[i+1:end]...;
							keep=keep)

end 

#===========================================================================#
#
# Given a relationship (item of type1) => (item of type2),
# Builds two dictionaries and returns two functions:
#				 - first gives the type-2 partner of a type1-item
#				 														(as the initial dict would do)
#				 - the second gives the type-1 partner of a type-2 item
#																		(as an inverse dict would do)
#
#
#---------------------------------------------------------------------------#


function FindPartners(pairs12;sortfirst=nothing)

	pairs12 = Dict(pairs12)

	isempty(pairs12) && return nothing,nothing

	pairs21 =	begin

							local K = vcat(keys(pairs12)...)
					
							isa(sortfirst,Bool) && sortfirst && sort!(K)

							isa(sortfirst,Function) && sort!(K, by=sortfirst)

							local V = unique(vcat(values(pairs12)...))

							Dict(v=>filter(k->in(v,vcat(pairs12[k])),K) for v in V)

						end



	return (key1 -> get(pairs12, key1, nothing), 
					key2 -> get(pairs21, key2, nothing))


#	return (key,label="") -> label==label1 ? pairs12[key] : pairs21[key]
end


#===========================================================================#
#
#	Quadrant of angle
#
#---------------------------------------------------------------------------#


function Quadrant(A::AbstractMatrix)

	size(A,2) == 2 && return Quadrant.(eachrow(A))

	size(A,1) == 2 && return Quadrant.(eachcol(A))

	error("Input not understood")

end

function Quadrant(xy::AbstractVector)

	length(xy) == 2 && return Quadrant(atan(xy[2],xy[1]))

	error("Input not understood")

end


function Quadrant(theta::Number)

	theta > pi && return Quadrant(theta-2pi)
	theta <-pi && return Quadrant(theta+2pi)

	for (i,(m,M)) in enumerate(([0,pi/2],[pi/2,pi],[-pi,-pi/2],[-pi/2,0]))
	
		m<=theta<=M && return i

	end


end


#===========================================================================#
#
# put "@time" if conditon==true
#
#---------------------------------------------------------------------------#

#macro timeif(cond,funwithargs)
#
#  @show funwithargs
#
#  cond && return @time eval(funwithargs)
#
#
#
#  return eval(funwithargs)
#
#end
#
#macro printif(cond,args...)
#
#  cond && println(join(args," "))
#
#end

#function timeif(cond)
#
#  cond==false && return f(expr,args...) = eval(expr)
#
#  return function f(expr,args...)
#
#    println(join(args," "))
#  
#    return @time eval(expr)
#
#  end
#
#
#end







#===========================================================================#
#
# convert from string to Symbol
#
#---------------------------------------------------------------------------#

function DictKey_Symbol(d)

  Dict([Symbol(k)=>v for (k,v) in pairs(d)])

end


#===========================================================================#
#
# Distribute list of parameters to several processes and print progress
#
#---------------------------------------------------------------------------#

function Distribute_Work(allparamcombs,do_work;arg_pos=1,kwargs0...)

  nr_scripts = min(length(allparamcombs), get_arg(1, arg_pos, Int64))

	start =  get_arg(1, arg_pos+1, Int64)

	idproc = gethostname()

	if start > nr_scripts

    println(string("\nI am ",idproc," and I am not doing any jobs (nr.jobs < nr.processes).\n"))
    return

  end

	stop = min(nr_scripts, get_arg(start, arg_pos+2, Int64))


	njobs = size(allparamcombs,1)

	which_ = sepLengths_cumulRanges(EqualDistributeBallsToBoxes(njobs, nr_scripts), start, stop)

	doparams = allparamcombs[which_] 
	


  
  doparams, which_, njobs = distribute_list(allparamcombs, nr_scripts, 
																						start, stop)
  
  println("\nI am $idproc and I am doing jobs $which_ out of $njobs.\n")


  function print_progress(ip,t1)

		println(string("\nI am $idproc and I completed $ip/",length(which_)," jobs (last one: ",which_[ip],"/$which_ in ", Int(round(Dates.value(Dates.now()-t1)/1000)),"s)."))

#    println(strout)

#    open("/home/pahomit/progress.txt","a") do fout
#
#      DlmF.writedlm(fout,[strout])
#
#    end

    return Dates.now()

  end



  for (ip,p) in enumerate(doparams)
  
    time1 = Dates.now()
  
    do_work(p;kwargs0...)

    print_progress(ip,time1)
  
  
  end


	return 


end


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function sepLengths_cumulRanges(L::AbstractVector{Int})::Vector{OrdinalRange{Int,Int}}
	
	boundaries = cumsum([0;L]) 
	
	return [boundaries[j]+1:boundaries[j+1] for j in axes(L,1)]

end 

function sepLengths_cumulRanges(L::AbstractVector{Int},
															 i::Int, j::Int=i)::OrdinalRange{Int,Int}

	i>j && error("'j'=$j must be at least 'i'=$i")

	range(1+sum(L[1:i-1]),step=1,length=sum(L[i:j]))

end 



function sepLengths_cumulRanges(As::Any)::Vector{OrdinalRange{Int,Int}}

	sepLengths_cumulRanges([size(A,1) for A in As])

end

#function sepLengths_cumulRanges(As...)::Vector{UnitRange}
#
#	sepLengths_cumulRanges(As)
#	
#end





#===========================================================================#
#
# Generate a list of n distinct random items from a list
#
#---------------------------------------------------------------------------#

function Random_Items(list, n=rand(axes(list,1)); distinct=true)

	if distinct
	
		n==length(list) && return Random.shuffle(list)

		n>length(list) && error("Cannot have so many distinct items")

		return list[Random.randperm(length(list))[1:n]]
	else 

		return [rand(list) for i=1:n]

	end 

#	function possible_extensions((L,N), c) 
#
#		[merge(c,Dict(length(c)+1=>i)) for i in Random.shuffle(axes(L,1))]
#
#	end 
#
#	function promising_candidate((L,N), c) 
#	
#		length(c)<=N && allunique(values(c))
#	
#	end 
#
#	accept_sol((L,N), c) = length(c)==N
#
#	function output(data, solutions, c)
#	
#		push!(solutions, c) 
#	
#		return false 
#	
#	end 
#
#	sol = Backtracking((list,n,distinct), 
#										 possible_extensions,
#										 promising_candidate,
#										 accept_sol,
#										 output
#										 )[1]
#
#	return [list[sol[i]] for i=1:n]



end


#===========================================================================#
#
# Write 
#
#---------------------------------------------------------------------------#






#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function is_dict_or_JLDAW(D)

	D isa AbstractDict && return true 

	all(in(propertynames(D)), [:keys, :values]) && return true

	return false 

end






#===========================================================================#
#
# Make directory tree if it doesn't exist
#
#---------------------------------------------------------------------------#

# function exsits already! 'mkpath'



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function logspace(start::Real, stop::Real, Len::Int64)

	exp.(range(log(start),log(stop),length=Len))

end

function uniqlogsp(start::Real, stop::Real, Len::Int64, tol::Real; Trunc=false)

	ntol,ftol = tolNF(tol)

	max_nr_steps = ftol>0 ? Int(floor(Float64(stop-start)/ftol)) : 100Len

	steps_step = max(1,Int(round((max_nr_steps-Len)/100.0)))

	for L in [Len:steps_step:max_nr_steps;max_nr_steps]

		S = Unique(logspace(start, stop, L), tol=tol) 

		if length(S)>=Len 
			
			s = S[Int.(round.(Rescale(1:Len,axes(S,1))))]

			return Trunc ? trunc.(s,digits=ntol) : s

		end 

	end 

	error("Interval too small or tolerance too high")


end 



#===========================================================================#
#
# Read variables from the arguments of the command
#
#---------------------------------------------------------------------------#

function get_arg(default,i,typ=String)

  if length(ARGS) >= i 

    typ!=String && return parse(typ,ARGS[i])

    return ARGS[i]

  end

  return default

end





#===========================================================================#
#
# "inverse map" --> applies the same function to one argument 
#
#---------------------------------------------------------------------------#




#===========================================================================#
#
# Executes a julia function in parallel
#
#---------------------------------------------------------------------------#




function execute_julia_function_(file_function;args=[],Nr_Procs=8,libs=[])


  lib_local  = "/media/tudor/Tudor/Work/scripts/julia_libraries/"
  lib_remote = "/net/horon/scratch/pahomit/apps/julia-1.2.0/lib/mylib/"

  lib_root = gethostname() == "tudor-HP" ? lib_local : lib_remote

		# make up some name for the script to be launched
  frun = string("auxrun_",join(file_function),rand(40000:50000),".jl")


  open(frun,"w") do f
    filename,funname = file_function

    text = string("include(\"",filename,".jl\")\n",funname,"()")
		# inclunde the file and run the target function

    write(f,text)
 #   println(text)
  end

  libs = "-L " *lib_root.*libs.*".jl"

  cmd = [["julia -p",Nr_Procs];libs;[frun];args]

  cmd = join(string.(cmd),' ')
		# give the input and output files and run in parallel

#  println(cmd)

  run(`sh -c $cmd`)

  rm(frun)

end



function execute_julia_function(args,nr_results,read_write_funs,file_function;Nr_Procs=8,libs=[])

  !isdir("./aux") && mkdir("./aux")

  write_args, read_result = read_write_funs



		# write the given arguments with the particular method
  f_args = write_args(args...)


		# make up some fileNames for the results
  f_result = "./aux/res_"*join(file_function).*string.(rand(20000:30000) .+ collect(1:nr_results) )

		# make a separate script and run the function in parallel 
  execute_julia_function_(file_function,args=[f_args;f_result],Nr_Procs=Nr_Procs,libs=libs)
		# read the results with the particular method
  return read_result(f_result)

end




#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#



function DictRandVals(params::T) where T<:AbstractDict

	T(k=>(typeof(v) <: AbstractVector) ? rand(v) : v 
												for (k,v) in pairs(params))

end

function DictRandVals(params::T) where T
	
	isList(T,AbstractDict) && return DictRandVals.(params)

	error("Type not understood: $T")

end

#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#



function DictFirstVals(params::T) where T<:AbstractDict

	T(k=>(typeof(v) <: AbstractVector) ? v[1] : v 
							for (k,v) in pairs(params))

end

function DictFirstVals(params::T) where T
	
	isList(T,AbstractDict) && return DictFirstVals.(params)

	error("Type not understood: $T")

end


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function AllValCombs(params::T;
										 constraint=nothing,sortby=nothing) where T<:AbstractDict

	vals(v) = (typeof(v) <: AbstractVector) ? v : [v]

	K = vcat(keys(params)...)

	dicts = [T(zip(K,vs)) 
					 for vs in Base.product([vals(params[k]) for k in K]...)][:]

  !isnothing(constraint) && filter!(constraint,dicts)

  !isnothing(sortby) && sort!(dicts,by=sortby)

  return dicts
end


function AllValCombs(params::AbstractVector{<:OrderedDict};kwargs...)

	collect(Base.product(AllValCombs.(params;kwargs...)...))[:]

end

#===========================================================================#
#
# Construct a NamedTuple from a list of keys and a dict/NT/array
#
#---------------------------------------------------------------------------#

function NT(data,Keys=nothing)


  if isa(data,AbstractDict) | isa(data,NamedTuple)

		Keys = Assign_Value(Keys,vcat(keys(data)...))

    data = [data[k] for k in Keys]
  end

#  data = try [data[k] for k in Keys]
#
#	 catch error
#
#           !any(isa.(error,[MethodError,ArgumentError])) && error("Please provide a list-type container (Array/...) or an associative collection (Dict/NamedTuple) with the right keys.")
#
#	   data
#
#	 end

  return NamedTuple{Tuple(Keys)}(data)

end


#===========================================================================#
#
# NamedTuple to Arrays (python compatibility)
#
#---------------------------------------------------------------------------#

function NT_toLists(NT)

  return [string(k) for k in keys(NT)], [NT[k] for k in keys(NT)]


end


#===========================================================================#
#
# For the common keys, the non-nothing values in 'nt' replace those in NT
#
#---------------------------------------------------------------------------#

function NT_ReplaceItems(NT,nt=NamedTuple{}())

  vals = Dict([k=>(k in keys(nt) ? nt[k] : NT[k]) for k in keys(NT)])

  K = Tuple([k for k in keys(NT) if !any(isnothing.(vals[k]))])

  return NamedTuple{K}([vals[k] for k in K])

end




function NT_KeepItems(NT,keep_keys=[])

  K = Tuple(intersect(keys(NT),keep_keys))

  return NamedTuple{K}([NT[k] for k in K])

end



function NT_AllCombs(params::NamedTuple;constraint=nothing,sortby=nothing)

  NTs = map(Base.product(values(params)...)) do vs

    return NamedTuple{keys(params)}(vs)
  end[:] 


  !isnothing(constraint) && filter!(constraint,NTs)

  !isnothing(sortby) && sort!(NTs,by=sortby)


  return NTs

    length(iter_last) == 0 && return NT_AllCombs_(params...)

end


#===========================================================================#
#
# Combine two lists of parameters (named tuples, dicts)
#
#---------------------------------------------------------------------------#

function Combine_NamedTuples(input_param,default_param)

  issubset(keys(input_param),keys(default_param)) || error("The parameters must be among "*join(map(string,keys(default_param)),", ")*".")

  return merge(default_param,input_param)

end




#===========================================================================#
#
# Pad a float left and right
#
#---------------------------------------------------------------------------#

function lrpad(x,N=[2,2];p='.')

	isa(N,Number) && return lrpad(x,[N,N];p=p)

  return rstrip(join([f(s[1:min(length(s),n)],n,'0') for (s,f,n) in zip(split(string(Float64(x)),'.'),[lpad,rpad],N)],p),p)

end


function nr2string(nr::T,digits=2) where T


	if T<:Union{Tuple,AbstractArray}# && any(t->typeof(nr[1])<:t,types)

		return map(n->nr2string(n,digits),nr) |> join
	
	end

	Union{Number,AbstractString} |> t-> T<:t || error("Please provide $t")



	isempty(digits) && return string(nr)

	digits isa Number && return nr2string(nr,[digits,digits])

  if nr isa AbstractString

    nr_ = tryparse(Float64,nr)
 
    return isnothing(nr_) ? nr : nr2string(nr_,digits)
 
  end


  left,right = split(string(round(Float64(abs(nr)),digits=digits[2])),'.')

  left = repeat("m",nr<0)*lpad(left,digits[1],'0')

  digits[2]==0 && return left

	right = rpad(right[1:min(end,digits[2])],digits[2],'0')

  return join([left,right],'p')

end

#===========================================================================#
#
#  Returns a path of n points which connects the inputed points
#
#---------------------------------------------------------------------------#

function PathConnect(points,n;end_point=true,bounds=[0,1],fdist=identity)

  size(points,1) == 1 && return points,[0]

  dist = LA.normalize(diff(points,dims=1) |>eachrow .|>LA.norm .|>fdist,1)

  n -= end_point


  ns_ = max.(1,Int64.(round.(dist*n)))

  while sum(ns_) > n 
    ns_[argmax(ns_)] -= 1
  end

  while sum(ns_) < n
    ns_[argmin(ns_)] += 1
  end
  

  path = vcat(map(axes(points,1)[1:end-1]) do i

    ep = end_point && i==axes(points,1)[end-1]

    t = range(0, 1, length = ns_[i] + 1 )[1:end-1+ep]

    return hcat(1 .- t , t)*points[i:i+1,:]

  end...)



  xticks = bounds[1] .+ cumsum([0;dist]).*diff(bounds)
#  xticks = cumsum([1:ns_])
 
  return path,xticks

end



#===========================================================================#
#
# return F(args) = f1(args) + f2(args) + ...
#
#---------------------------------------------------------------------------#

function Sum_functions(functions...)::Function
 

#  functions = [f for f in functions if isa(f,Function)]

#  return (args...) -> mapreduce(F -> F(args...),+, Functions)
  return (args...) -> reduce(+,map(f -> f(args...), functions))

end


#===========================================================================#
#
# Make array out of list of lists, like np.array(...)
#
#---------------------------------------------------------------------------#

function subarray(A,js=[])

  (isa(A,Number) | (length(js)==0)) && return A

  return subarray(A[js[1],axes(A)[2:ndims(A)]...],js[2:length(js)])

end


function get_shape(A)

  shape = Int64[]

  while true
  
    isa(A,Number) && return Tuple(shape)

    push!(shape,size(A,1))

    A = subarray(A,[1])
   
  end

end

	# much faster than my function
#function nparray(L)::Array
#
#  np = PyCall.pyimport("numpy")
#
#  return np.array(L)
#
#end

#function ListToArray(L,dtype=Float64)
#
#  A = zeros(dtype,get_shape(L))
#
#  for i in CartesianIndices(A)
#
#    A[i] = subarray(L,Tuple(i))
#
#  end  
#
#  return A
#
#end

#===========================================================================#
#
# vectors of integers, each component is a range
#
#---------------------------------------------------------------------------#

function vectors_of_integers(D::Int64, stop, start=-stop; dim=1, sortby=nothing)

	dim2 = [2,1][dim]


	boundaries = zipmap([Assign_Value(start,-stop), stop]) do x
		
		isa(x,Int) && return  fill(x,D)
	
		isList(x,Int) || error("Type ",typeof(x)," not supported")

		length(x)>=D && return vcat(x...)[1:D]

		length(x)==1 && return fill(x[1], D)

		error("Wrong dimensions")
	
	end 


	Ls = [1;[f-i+1 for (i,f) in boundaries];1]


	out = similar(Array{Int64}, [prod(Ls), D][[dim,dim2]]...)

	
	for (i,(a,b)) in enumerate(boundaries)

		setindex!(out,
							repeat(a:b, inner=prod(Ls[i+2:end]), outer=prod(Ls[1:i])),
							[:,i][dim], [:,i][dim2]
							)


	end

	!isa(sortby, Function) && return out

	return sortslices(out, dims=dim, by=sortby)





end





#===========================================================================#
#
# Check if we have the same number or the same array
#
#---------------------------------------------------------------------------#


function fSame(atol)

  function same(a::Number,b::Number=0.0,atol=atol)::Bool

    return isapprox(a,b,atol=atol)

  end



  function same(a::AbstractArray,b::AbstractArray,atol=atol)::Bool

    return isapprox(a,b,atol=atol)

  end



  function same(a::AbstractArray,b::Number=0.0,atol=atol)::Bool

    return isapprox(LA.norm(a),b,atol=atol)
  
  end



  function same(a::Number,b::AbstractArray,atol=atol)::Bool

    return isapprox(LA.norm(b),a,atol=atol)
  
  end


  return same

end



#function Same(a,b,order::Real=6)::Bool
#
#  return LA.norm(a.-b) < 1.0/10.0^order
#
#end




#===========================================================================#
#
# Assigns value if variable is nothing
#
#---------------------------------------------------------------------------#


function Assign_Value(input_value, std_value)

  isnothing(input_value) && !isnothing(std_value) && return std_value
  
  return input_value

end


function Assign_Value(input_value, get_std_value, args...; kwargs...)

	!isnothing(input_value) && return input_value
	
	!isnothing(get_std_value) && return get_std_value(args...; kwargs...)

  return nothing

end



#############################################################################

end
