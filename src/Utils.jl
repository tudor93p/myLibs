module Utils
#############################################################################



import ..LA 


import Dates, Combinatorics


using Distributed 

using OrderedCollections:OrderedDict
import Random


const List = Union{AbstractVector, AbstractSet, Tuple, Base.Generator}




#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

#applyMethod(f::Function, methods::List,  




#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function Backtracking(data,
											possible_extensions::Function, 
											promising_candidate::Function, 
											accept_sol::Function,
											root::Function,
											solutions=[],
											)::Vector

	Backtracking(data,
							 possible_extensions, 
							 promising_candidate, 
							 accept_sol, 
							 (data, solutions, candidate) -> push!(solutions, candidate),
							 root,
							 solutions
							)
end 

function Backtracking(data,
											possible_extensions::Function, 
											promising_candidate::Function, 
											accept_sol::Function, 
											output::Function,
											root::Function,
											solutions=[],
											)::Vector


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


	backtrack!(data, solutions, root(data))

	return solutions 

end 








#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#






function DistributeBallsToBoxes(balls::Int, boxes::Int, 
																expand::Function,
															 )::Vector{Vector{Int}}
	
	# 'expand' is the inverse of 'effval' 
	# effval(xi)==0 for xi in expand(x) 
	#	effval: many to one, expand: one to many 

	boxes==0 && return [Int[]]

	return Utils.flatmap(DistributeBallsToBoxes_(UInt(balls), boxes)) do X

		iter = Base.product((unique(vcat(expand(x)...)) for x in X)...)

		return collect(Iterators.partition(Iterators.flatten(iter),boxes))

	end 

end 






function DistributeBallsToBoxes(balls::Int, boxes::Int)::Vector{Vector{Int}}

	DistributeBallsToBoxes_(UInt(balls), boxes)

end 






function DistributeBallsToBoxes_(balls::UInt, boxes::Int
															 )::Vector{Vector{Int}}

	boxes==0 && return [Int[]]

	return map(Combinatorics.combinations(1:(balls+boxes-1), boxes-1)) do d

		convert(Vector{Int}, diff(vcat(0, d, balls+boxes)) .- 1)

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


function PropDistributeBallsToBoxes(balls::Int, sizes::AbstractVector{<:Number})::Vector{Int} 

	float_ns = balls*LA.normalize(sizes, 1) # proportional distribution 

	int_ns = Int64.(round.(float_ns))		# rough distribution 


	while sum(int_ns) > balls 

		nonzero_inds = findall(>(0), int_ns)

		isempty(nonzero_inds) && return int_ns 

		jmax = argmax(int_ns[nonzero_inds] - float_ns[nonzero_inds])

		int_ns[nonzero_inds[jmax]] -= 1 

  end

  while sum(int_ns) < balls 

		int_ns[argmin(int_ns - float_ns)] += 1

  end

	return int_ns 

end 




function EqualDistributeBallsToBoxes_cumulRanges(args...)

	d = EqualDistributeBallsToBoxes(args[1:2]...)

	return sepLengths_cumulRanges(d, args[3:end]...)

end 


function PropDistributeBallsToBoxes_cumulRanges(args...)

	d = PropDistributeBallsToBoxes(args[1:2]...)

	return sepLengths_cumulRanges(d, args[3:end]...)

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



function Unique!(V::AbstractArray; sorted=false, kwargs...)

	i = Unique(V; kwargs..., sorted=false, inds=:first)[2]

	deleteat!(V, filter(!in(i), axes(V,1)))

	sorted && sort!(V)

end



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#



tolNF(tol::Int)::Tuple{Int,Float64}  = (tol,10.0^(-tol))

function tolNF(tol::Float64)::Tuple{Int,Float64}

	tol==0.0 ? 50 : Int(ceil(-log10(tol))), tol

end 




#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


#function valofkey_dictorJLDAW(data::Union{AbstractDict,NamedTuple}, k)
#
#	get(data, k, nothing)
#
#end 
#
#
#
#function valofkey_dictorJLDAW(JLDAW_object, k)
#
#	for (i,ki) in enumerate(JLDAW_object.keys)
#
#		k==ki && return JLDAW_object.values[i]
#		
#	end 
#	
#	return nothing 
#
#end  





#function gk_gv_vk_dictorJLDAW(data::Union{AbstractDict,
#																					NamedTuple})::NTuple{3,Function}
#
#	(keys, values, valofkey_dictorJLDAW)
#
#end  
#
#
#
#
#function gk_gv_vk_dictorJLDAW(data)::NTuple{3,Function}
#	
#	(d->d.keys, d->d.values, valofkey_dictorJLDAW) 
#
#end  
#



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#




function Unique(V::AbstractArray{T};
								dim=1,
								tol=1e-8, inds=nothing,
								sorted=false,
								check_type=true) where T

	if !check_type || is_exact(T) || all(is_exact ∘ typeof, V)

		U = unique(V,dims=dim) |> function (u)

				sorted isa Bool && !sorted && return u 

				kw = sorted isa Bool ? () : Dict(:by=>sorted)

				ndims(V)==1 && return sort(u; kw...)
	
				return sort(u; dims=dim, kw...)	

			end 


		isnothing(inds) && return U

		i_f = findfirst(isequal(inds), [:first,:all,:last])

		isnothing(i_f) && error("Cannot get indices with '$inds'")

		f = [findfirst, findall, findlast][i_f]

		iter = collect(eachslice(V, dims=dim))

		return U, [f(isequal(u), iter) for u in eachslice(U,dims=dim)]

	end 


	 
	function get_newV()

		ntol,ftol = tolNF(tol)

		if ftol>0 && (is_float(T) || any(is_float ∘ typeof, V))
			
			return trunc.(V, digits=ntol)

		end 

		return V

	end


	I = Unique(get_newV(); tol=tol, inds=Assign_Value(inds, :first), dim=dim,
						 						 sorted=sorted, check_type=false)[2]

	isnothing(inds) && return V[I]

	return V[first.(I)], I


end 

Unique(V::Tuple; kwargs...) = Unique(collect(V); kwargs...)


function EnumUnique(V; kwargs...)
	
	zip(Unique(V; kwargs..., inds=:all)...)

end 

#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#



function flatmap(f::Function, itr)#; kwargs...)

#	mapreduce(f, vcat, itr, init=[])

	vcat(map(f,itr)...)

end

function flatmapif(f::Function, pred, itr)

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

	filter(pred, map(f,itr))

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



invmap(arg, fs::Vararg{<:Function}; kwargs...) = invmap(arg, fs; kwargs...)

function invmap(args, fs::List; kwargs...)

	map(f->f(args...; kwargs...), fs)

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

function IdentifySectors(f::Function, iter::AbstractVector, 
												 args...; kwargs...)::Vector{UnitRange{Int}}

	IdentifySectors(map(f, iter), args...; kwargs...)

end 

function IdentifySectors(list::AbstractVector{T}; kwargs...
												 )::Vector{UnitRange{Int}} where T<:Union{
																							<:Int, <:AbstractSet, <:Bool,
																							<:Tuple{Vararg{<:Int}}, 
																							<:Tuple{Vararg{<:Bool}},
																							<:AbstractVector{<:Int}, 
																							<:AbstractVector{<:Bool},
																	}
	IdentifySectors_(==, list)

end 


function IdentifySectors(list::AbstractVector{T}; tol=1e-8, kwargs...
												 )::Vector{UnitRange{Int}} where T<:Union{
														<:ComplexF64, <:Float64,
														<:AbstractVector{<:Union{<:ComplexF64,<:Float64}},
																	}

	ntol,ftol = tolNF(tol)
	
	return IdentifySectors_(list) do a,b

		isapprox(a,b, atol=ftol)

	end 

end 





function IdentifySectors_(same::Function, list::AbstractVector, 
												 sectors::Vector{UnitRange{Int}}=UnitRange{Int}[], 
												 start::Int=0
												 )::Vector{UnitRange{Int}} 

	update(j) = vcat(sectors, UnitRange{Int}[1+start:j+start])


	for i in axes(list,1)

		if same(list[i], list[1])		# still within the sector
		
			i==length(list) && return update(i) # sector ended abruptly

		else  # sector ended already at i-1

			return IdentifySectors_(same, list[i:end], update(i-1), start+i-1)

		end
	end

end




#===========================================================================#
#
#	Identify Ranges -- group together inds which form ranges 
#
#---------------------------------------------------------------------------#


function IdentifyRanges(list::AbstractVector{<:Int}
												)::Vector{Union{Int,OrdinalRange{Int,Int}}}

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

	sector(i) = Assign_Value(findfirst(s->in(i,s),sectors),0)

	return flatmap(IdentifySectors(sector, axes(list,1))) do S

		!in(S,sectors) && return list[S]

		(a,b) = extrema(S)
		
		step = list[a+1]-list[a]

		step==1 && return [UnitRange(list[a], list[b])]
		
		return [StepRange(list[a], step, list[b])]

	end


end




#===========================================================================#
#
# check if arg is a vector or a tuple 
#
#---------------------------------------------------------------------------#




function isTuple(arg::T, Ti=Any)::Bool where T

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



function isList(arg::T, Ti=Any)::Bool where T

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


	return isTuple(isa(arg,Base.Generator) ? Tuple(arg) : arg, Ti)

end


function isList(;T=Any)::Function 
	
	islist(arg)::Bool = isList(arg, T)
	
end 



function isListLen(X, L::Int)::Bool 

	isList(X) && length(X)==L

end 


function isListLen(X, T, L::Int)::Bool 

	isList(X, T) && length(X)==L

end 



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function flat(list_of_lists...; keep=nothing, flatten=[List])::Vector 

	i = findfirst(list_of_lists) do x 

				any(flatten) do T

					x isa T

				end 

			end 

	if isnothing(i) 
	
		keep isa Function && return filter(keep, vcat(list_of_lists...))

#		isnothing(keep) &&  
		return vcat(list_of_lists...)

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

function Quadrant(A::AbstractMatrix, C::AbstractMatrix; dim::Int)::Int

	Quadrant(A.-C; dim=dim)

end

function Quadrant(A::AbstractMatrix; dim::Int)::Int

	Quadrant.(eachslice(A, dims=dim))

end

function Quadrant(xy::AbstractVector, C::AbstractVector)::Int

	Quadrant(xy-C)

end 

function Quadrant(xy::AbstractVector)::Int

	length(xy) == 2 && return Quadrant(atan(xy[2],xy[1]))

	error("Input not understood")

end


function Quadrant(theta::Number)::Int

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

function Distribute_Work(allparamcombs::AbstractVector,
												 do_work::Function; 
												 vararg::Bool=false, arg_pos::Int=1, kwargs0...)

  nr_scripts = min(length(allparamcombs), get_arg(1, arg_pos, Int64))

	start =  get_arg(1, arg_pos+1, Int64)

	idproc = gethostname()

	if start > nr_scripts

    println(string("\nI am $idproc and I am not doing any jobs (nr.jobs < nr.processes).\n"))
    return

  end

	stop = min(nr_scripts, get_arg(start, arg_pos+2, Int64))


	njobs = size(allparamcombs,1)

	which_ = EqualDistributeBallsToBoxes_cumulRanges(njobs, nr_scripts, start, stop)

	doparams = allparamcombs[which_] 
	

  
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
 
		vararg ? do_work(p...; kwargs0...) : do_work(p; kwargs0...)

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

	return range(1+sum(L[1:i-1]),step=1,length=sum(L[i:j]))

end 

#sepLengths_cumulRanges(L)[i:j+1]

function sepLengths_cumulRanges(As::Any)::Vector{OrdinalRange{Int,Int}}

	sepLengths_cumulRanges([size(A,1) for A in As])

end






#===========================================================================#
#
# Generate a list of n distinct random items from a list
#
#---------------------------------------------------------------------------#

function Random_Items(list::List, n=rand(axes(list,1)); distinct=true)::List

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

#is_dict_or_JLDAW(D::Union{NamedTuple,AbstractDict})::Bool = true 
#
#function is_dict_or_JLDAW(D)::Bool
#
#	all(in(propertynames(D)), [:keys, :values]) && return true
#
#	return false 
#
#end






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

function Rescale(A::AbstractArray{<:Number,N}, mM0, mM1=A)::AbstractArray{<:Number,N} where N 

  m0,M0 = extrema(mM0)

	length(A)==1 && return fill(m0, size(A))

	m,M = extrema(mM1)
	
	return (A .- m)*(M0-m0)/(M-m) .+ m0

end




#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


RescaleInds(N::Int, args...)::Vector{Int} = RescaleInds(1:N, args...)

function RescaleInds(A::AbstractArray{<:Number}, d::Int,
											large_arg1::Union{Int,AbstractArray{<:Number}},
											large_args...)::Vector{Int}
	
	RescaleInds(axes(A, d), large_arg1, large_args...)

end 


function RescaleInds(A::AbstractVector{Int},
											large_number::Int)::Vector{Int} 

	RescaleInds(A, 1:large_number)

end 

function RescaleInds(A::AbstractVector{Int},
											large_array::AbstractArray{<:Number}, d::Int
											)::Vector{Int} 

	RescaleInds(A, axes(large_array, d))

end 


function RescaleInds(small_range::Ts, large_range::Tl
											)::Vector{Int} where {T<:AbstractVector{Int},
																						Ts<:T, Tl<:T}

	@assert issubset(small_range, UnitRange(extrema(large_range)...))

	out = Int.(round.(Rescale(small_range, large_range)))

	@assert length(unique(out))==length(small_range)

	return out 

end  









#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function linspace(start::Real, stop::Real, Len::Int64)::Vector{Float64}

	collect(range(start, stop, length=Len))

end 

function logspace(start::Real, stop::Real, Len::Int64)::Vector{Float64}

	exp.(linspace(log(start), log(stop), Len))

end






function uniqsp(getsp::Function)::Function 
	
	function uniqsp_(start::Real, stop::Real, Len::Int64, tol::Real; Trunc=false)::Vector{Float64}

		ntol,ftol = tolNF(tol)
	
		max_nr_steps = ftol>0 ? Int(floor(Float64(stop-start)/ftol)) : 100Len
	
		steps_step = max(1,Int(round((max_nr_steps-Len)/100.0)))
	
		for L in [Len:steps_step:max_nr_steps;max_nr_steps]
	
			S = Unique(getsp(start, stop, L), tol=tol) 
	
			if length(S)>=Len 
				
				s = S[RescaleInds(Len, S, 1)]

				return Trunc ? trunc.(s,digits=ntol) : s
	
			end 
	
		end 
	
		error("Interval too small or tolerance too high")
	end 

end  



uniqlogsp = uniqsp(logspace) 

uniqlinsp = uniqsp(linspace)

#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#





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


dict_constr(::OrderedDict) = OrderedDict
dict_constr(::AbstractDict) = Dict 
dict_constr(::NamedTuple) = NamedTuple 


function dict_like(source::Union{<:AbstractDict, <:NamedTuple}, dest::List)
	
	dict_constr(source)(dest)
	
end 
	
function dict_like(source::Union{<:AbstractDict, <:NamedTuple},
									 dest::Union{<:AbstractDict, <:NamedTuple}=source)

	dict_like(source, (k=>v for (k,v) in pairs(dest)))

end 


function dict_like(source::Union{<:AbstractDict, <:NamedTuple}, dest::Pair)

	dict_like(source, [dest])

end 


function dict_like(source::Union{<:AbstractDict, <:NamedTuple}, 
									 dest::Base.Iterators.Zip)

	dict_like(source, (k=>v for (k,v) in dest))

end 



function adapt_merge(D0::Union{<:AbstractDict, <:NamedTuple}, dicts...)

	merge(D0, (dict_like(D0, d) for d in dicts)...) 

end 



getkeys(a::AbstractVector)::Vector = a
getkeys(a::Union{<:AbstractDict, <:NamedTuple})::Vector = collect(keys(a))
getkeys(a::Pair)::Vector = [p.first]
getkeys(a::Union{<:Int,<:Tuple,<:Symbol,<:AbstractString})::Vector = [a]


function dict_diff(D::Union{<:AbstractDict, <:NamedTuple},
									 args...)
	
	K = setdiff(keys(D), vcat(getkeys.(args)...))

	return dict_like(D, (k=>D[k] for k in K))

end 


function dict_keepkeys(D::Union{<:AbstractDict, <:NamedTuple},
											 args...)
	
	K = intersect(keys(D), vcat(getkeys.(args)...))
	
	return dict_like(D, (k=>D[k] for k in K))

end 




#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#



function pickrand(params)

	pick(v::AbstractVector) = rand(v)
	pick(v) = v 

	return (k=>pick(v) for (k,v) in params)

end 








function DictRandVals(params::Union{<:AbstractDict, <:NamedTuple}
										 )::Union{<:AbstractDict, <:NamedTuple}

	dict_constr(params)(pickrand(pairs(params)))

end




function DictRandVals(params::List)

	#isList(params, Union{<:AbstractDict,<:NamedTuple}) && 
	return DictRandVals.(params)

	#error("Type not understood: ",typeof(params))

end

#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#



function DictFirstVals(params::Union{<:AbstractDict, <:NamedTuple}
										 )::Union{<:AbstractDict, <:NamedTuple}

	dict_constr(params)(k=>isa(v,AbstractVector) ? v[1] : v 
							for (k,v) in pairs(params))

end

function DictFirstVals(params::List)
	
	isList(params, Union{<:AbstractDict, <:NamedTuple}) && return DictFirstVals.(params)

	error("Type not understood: ",typeof(params))

end


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function AllValCombs(params::Union{<:NamedTuple, <:AbstractDict};
										 constraint=nothing, sortby=nothing
										 )::Vector{<:Union{<:NamedTuple, <:AbstractDict}}

	vals(v::AbstractVector) = v 
	vals(v) = [v]
	
	K = collect(keys(params))

	constr = dict_constr(params)

	dicts = [constr(zip(K,vs)) 
					 for vs in Base.product((vals(params[k]) for k in K)...)][:]

  isnothing(constraint) || filter!(constraint, dicts)

  isnothing(sortby) || sort!(dicts, by=sortby)

  return dicts 

end


function AllValCombs(params::List; kwargs...)

	collect(Base.product(AllValCombs.(params; kwargs...)...))[:]

end

#===========================================================================#
#
# Construct a NamedTuple from a list of keys and a dict/NT/array
#
#---------------------------------------------------------------------------#

function NT(data::Union{<:AbstractDict,<:OrderedDict,<:NamedTuple},
						Keys=collect(keys(data)))::NamedTuple

  NamedTuple{Tuple(Keys)}([data[k] for k in Keys])

end 


#function NT(data::List, args...)::Vector{NamedTuple}
#	
#	[NT(d) for d in data]
#
#end 


#===========================================================================#
#
# NamedTuple to Arrays (python compatibility)
#
#---------------------------------------------------------------------------#

function NT_toLists(NT)

	[string(k) for k in keys(NT)], [NT[k] for k in keys(NT)]

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

function lrpad(x::Real, N::Number; kwargs...)::String

	lrpad(x, [N,N]; kwargs...)

end 

function lrpad(x::Real, N::AbstractVector{Int}=[2,2]; p='.')::String 

  rstrip(join([f(s[1:min(length(s),n)],n,'0') for (s,f,n) in zip(split(string(Float64(x)),'.'),[lpad,rpad],N)],p),p)

end



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#



function nr2string(nr::List, args...)::Vector{String}

	[nr2string(n, args...) for n in nr]

end 




function nr2string(nr::Real, digits::Number=2)::String 

	nr2string(nr, [digits] )

end 

function nr2string(nr::Real, digits::Tuple)::String 

	nr2string(nr, convert(Vector{Int}, collect(digits)))

end 


function nr2string(nr::AbstractString, args...)::String

	nr_ = tryparse(Float64, nr)

	return isnothing(nr_) ? nr : nr2string(nr_, args...)
	
end 



function nr2string(nr::Union{Real,AbstractString}, digits::AbstractVector{Int})::String 

	isempty(digits) && return string(nr) 

	length(digits)==1 && return nr2string(nr, [digits; 0]) 


  left,right = split(string(round(Float64(abs(nr)),digits=digits[2])),'.')

  left = repeat("m",nr<0)*lpad(left,digits[1],'0')

  digits[2]==0 && return left

	right = rpad(right[1:min(end,digits[2])],digits[2],'0')

	return string(left,'p',right)

end


#===========================================================================#
#
# Recast a vector as a single-column (dim=2) or single-row matrix (dim=1) 
# 	works with any-dimensional arrays
#
#---------------------------------------------------------------------------#

function new_axis(nr_dims::Int, dim::Int)::Int 

	dim==1 && return 1

	dim==2 && return nr_dims+1

	error()

end 

function new_axis(x::AbstractArray{T,N}, dim::Int)::Int where {T,N}

	new_axis(N, dim)

end 






function Add1Dim(A::AbstractArray{T,N}, dim::Int, new_dim::Int
												 )::AbstractArray{T,N+1} where T<:Number where N

	reshape(A, insert!([size(A)...], new_dim, 1)...) 

end 


function Add1Dim(A::AbstractArray{T,N}, dim::Int
												 )::AbstractArray{T,N+1} where T<:Number where N

	Add1Dim(A, dim, new_axis(A, dim))

end

function Add1Dim(dim::Int)::Function

	function a1d(A::AbstractArray{T,N}
							 )::AbstractArray{T,N+1} where T<:Number where N

		Add1Dim(A, dim)

	end 

end 


function Add1Dim(nr_dims::Int, dim::Int)::Function

	a = new_axis(nr_dims, dim) 

	function a1d(A::AbstractArray{T,N}
							 )::AbstractArray{T,N+1} where T<:Number where N

		@assert N==nr_dims

		return Add1Dim(A, dim, a)

	end 

end 






function VecAsMat(V::AbstractArray{T}, dim::Int
								 )::AbstractMatrix{T} where T<:Number

	count(size(V).>1)>1 && error("Wrong shape, vector not understood")

	return Add1Dim(view(V, :), dim)

end 


function VecAsMat(dim::Int)::Function 

	a = new_axis(1, dim) 

	return function vam(V::AbstractArray)::AbstractMatrix

		Add1Dim(view(V,:), dim, a) 

	end 

end 





function sel(dim::Int)::Function 

	function out(A::AbstractArray, inds)

		selectdim(A, dim, inds) 

	end 

end 


function sel(dim::Int, inds)::Function 

	function out(A::AbstractArray)  

		selectdim(A, dim, inds) 

	end 

end 


function sel(A::AbstractArray, dim::Int)::Function

	function out(inds)

		selectdim(A, dim, inds)

	end 

end 


#===========================================================================#
#
#  Returns a path of n points which connects the inputed points
#
#---------------------------------------------------------------------------#



function PathConnect(points::AbstractVector, n::Int;
										 kwargs...)::Tuple{Vector,Vector{Float64}} 

	PathConnect(points, n, LA.norm.(diff(points)); kwargs...)

end 

function PathConnect(points::AbstractVector, n::Int, dist::AbstractVector;
										 kwargs...)::Tuple{Vector,Vector{Float64}}

	reshape.(PathConnect(VecAsMat(points, 1), n, dist; dim=2, kwargs...), :)

end 


#										 fdist::Function=identity,

function PathConnect(points::AbstractMatrix, n::Int; 
										 dim::Int=1, kwargs...)::Tuple{Matrix,Vector{Float64}}

	PathConnect(points, n,
							LA.norm.(eachslice(diff(points,dims=dim),dims=dim));
							dim=dim, kwargs...)

end 




function PathConnect(points::AbstractMatrix, n::Int, dist::AbstractVector;
										 end_point::Bool=true, 
										 bounds::AbstractVector=[0,1],
										 dim=1,
										 )

	nr_intervals = size(points, dim) - 1

	nr_intervals > 0 || return points, [] 

	nr_intervals==length(dist) || error("There should be $nr_intervals distance(s), you provided ",length(dist))

	xticks = Rescale(cumsum([0;dist]), bounds) 


	if n <= nr_intervals + end_point 
	
		return selectdim(points, dim, 1:n), xticks[1:n]

	end 
	

	ns = PropDistributeBallsToBoxes(n-end_point, dist)  


	for i in findall(ns .== 0)[1:min(end,count(ns.>0))]

		ns[argmax(ns)] -= 1 # "steal" one unit 

		ns[i] += 1

  end
 

	dim2 = [2,1][dim] 

	path = mapreduce([vcat,hcat][dim], 1:nr_intervals) do i

    ep = end_point && i==nr_intervals

		t = range(0, 1, length = ns[i] + 1)[1:end-!ep] |> VecAsMat(dim2)
	
		return CombsOfVecs(selectdim(points, dim, i:i+1), 
											 cat(1 .- t , t, dims=dim2), dim=dim)

	end


  return path, xticks

end



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


#merge(D1, dict_constr(D1)(pairs(D2)))



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function catNewAx(dim::Int, arrays::Vararg{<:AbstractArray})::Array{<:Number}

	Ns = ndims.(arrays)

	N = minimum(Ns)


	@assert issubset(Ns, [N, N+1])
	
	d = new_axis(N, dim)  


	shapes = map(zip(arrays,Ns)) do (a,n)

		inds = collect(1:n)

		n==N+1 && setdiff!(inds, d)

		return size(a)[inds]

	end 

	@assert length(unique(shapes))==1 


	function prep(a::AbstractArray{t,n}
								)::AbstractArray{t} where n where t<:Number
	
		n==N ? Add1Dim(a, dim, d) : a

	end 


	return cat((prep(a) for a in arrays)..., dims=d)

end  



function VecsToMat(vecs::Vararg{AbstractArray}; dim::Int)::Matrix 

	cat((VecAsMat(v, dim) for v in vecs)..., dims=dim)

end 





#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#



function RecursiveMerge_(;kwargs...)::Function  

	(args...) -> RecursiveMerge_(args...; kwargs...)

end 


function RecursiveMerge_(x::Union{<:Number,<:AbstractVector{<:Number}},
												 y::Number;
												 kwargs...)::Vector{<:Number}

	vcat(x,y)

end 


function RecursiveMerge_(u::AbstractArray{<:Number,N},
												 v::AbstractArray{<:Number,N}; 
												 dim::Int, kwargs...
												 )::Array{<:Number, N+1} where N 

	@assert N in [1,2]

	return catNewAx(dim, u, v)	

end 


function RecursiveMerge_(u::AbstractArray{<:Number,N},
												 v::AbstractArray{<:Number,M}; 
												 dim::Int, kwargs...
												 )::Array{<:Number, N} where {N,M}

	@assert N==M+1 
	
	N>3 && @warn "Strange array"

	return catNewAx(dim, u, v)	

end 






function RecursiveMerge_(A::AbstractDict{Ka,Va}, B::AbstractDict{Kb,Vb}; 
												 kwargs...)::AbstractDict{promote_type(Ka,Kb), Any} where {Ka,Kb,Va,Vb}
	
	RecursiveMerge(A, B; kwargs...)

end 


function RecursiveMerge(arg1::AbstractDict, args...; 
												kwargs...)::Dict{<:Any, Any} 

	# mergewith does not preserve OrderedDict! It turns it into Dict 
	# hence empty dict needed 
	
	merge(RecursiveMerge_(;kwargs...), empty(arg1, Any), arg1, args...)

end  


function RecursiveMerge(; kwargs...)::Function 

	(args...) -> RecursiveMerge(args...; kwargs...)

end  


function RecursiveMerge(arg::AbstractDict{K,V}; kwargs...)::AbstractDict{K,V} where {K,V}

	arg 

end  

function RecursiveMerge(arg::List; kwargs...)::AbstractDict{<:Any, Any}

	isList(arg, AbstractDict) || error("Argument not understood")

	return RecursiveMerge(arg...; kwargs...)

end 

	

function RecursiveMerge(f::Function, iter::List; parallel=false, kwargs...)::AbstractDict{<:Any, Any}

	RecursiveMerge((parallel ? pmap : map)(f, iter); kwargs...)

end 




#===========================================================================#
#
# return F(args) = f1(args) + f2(args) + ...
#
#---------------------------------------------------------------------------#

function Sum_functions(functions::Vararg{Function})::Function
 
#  return (args...) -> mapreduce(F -> F(args...),+, Functions) 

	(args...) -> reduce(+, invmap(args, functions))

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
#
#
#---------------------------------------------------------------------------#

function CombsOfVecs(vecs::AbstractMatrix{Tv}, 
										 coeff::AbstractMatrix{Tc}; dim::Int
										 )::Matrix{promote_type(Tv,Tc)} where {Tv,Tc}

	dim==2 && return vecs*coeff # Vectors are on columns 

	dim==1 && return coeff*vecs # the vectors are on lines 

	error()

end


function CombsOfVecs10(A::AbstractMatrix{<:Number}, 
											 stopstart...; dim::Int)::Matrix{<:Number}

	CombsOfVecs(A,
							vectors_of_integers(size(A,dim), stopstart...; dim=dim),
							dim=dim)

end 




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
#
#
#---------------------------------------------------------------------------#

function getprop(prop::Union{List,Symbol}, std_value=nothing)::Function

	M -> getprop(M, prop, std_value)

end 


function getprop(M, prop::Symbol, std_value=nothing)

	allnames = isa(M,Module) ? names(M, all=true, imported=true) : propertynames(M)

	return in(prop, allnames) ? getproperty(M, prop) : std_value 

end 


function getprop(M, prop::List, std_value=nothing)

	std = isList(std_value) ? i->std_value[i] : i->std_value

	return [getprop(M, p, std(i)) for (i,p) in enumerate(prop)]
						
end 



	

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


#===========================================================================#
#
# Split tasks on workstations according to difficulty (cost)
#
#---------------------------------------------------------------------------#


#############################################################################
module SplitWorkstations
#############################################################################


import ..Assign_Value, ..OrderedDict, ..IdentifySectors, ..Backtracking
import Combinatorics

#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#



function cost(M::AbstractVector{<:Real}, 
							possible_machines::AbstractVector{<:AbstractDict},
							inds::AbstractVector{Int},
							machines::AbstractVector{<:AbstractString})::Float64

	sum(M[i]/sum(possible_machines[i][m] for m in machines) for i in inds)

end  


function cost_constraint(C::Union{<:AbstractVector{<:Real},
																	<:Tuple{Vararg{<:Real}}};
												 at_most_from_prev::Real=1.25,
												 kwargs...)::Bool 

	# rough decreasing order  C[i+1] < factor*C[i]  , factor close to 1

	length(C)<=1 && return true  

	return all(diff(C)./C[1:end-1] .+ 1 .< at_most_from_prev)
	
end  

cost_constraint(C::Vararg{<:Real}; kwargs...)::Bool = cost_constraint(vcat(C...); kwargs...)




function cost_constraint2(C::Union{<:AbstractVector{<:Real},
																	<:Tuple{Vararg{<:Real}}};
													at_least_from_prev::Real=0.6,
													kwargs...)::Bool 
	

	length(C)<=1 && return true 

	return all(diff(C)./C[1:end-1] .+ 1 .>  at_least_from_prev)


end 

cost_constraint2(C::Vararg{<:Real}; kwargs...)::Bool = cost_constraint2(vcat(C...); kwargs...)




struct IterSplitWS 

	possible_machines::Vector{OrderedDict{String,Int}}

	sectors::Vector{Vector{Int}}	

	M::Vector{Float64}

	m::Vector{String}  

	min_ms::Int

	start::Int

	i::Int
	
	cconstraint::Function 

end  

function IterSplitWS(data::Tuple{<:Any,<:Any,<:Any},
						 m::AbstractVector{<:AbstractString}, 
						 min_ms::Int=length(m), 
						 start::Int=1, 
						 i::Int=1;
						 kwargs...)

	IterSplitWS(data..., m, min_ms, start, i, C->true)

end 

function IterSplitWS(data::Tuple{<:Any,<:Any,<:Any},
						 m::AbstractVector{<:AbstractString}, 
						 min_ms::Int, 
						 start::Int, 
						 i::Int,
						 cconstraint::Real;
						 kwargs...)


	IterSplitWS(data...,
			m, min_ms, start, i,
			C->cost_constraint(cconstraint, C; kwargs...),
			)

end 




#===========================================================================#
#
# Iterate to find possible extensions 
#
#---------------------------------------------------------------------------#




function next_m(possible_machines, inds,
								all_m::AbstractVector{<:AbstractString}, 
								min_m::Int, outside_m::Int
								)::Tuple{<:Any,<:Any}

	A = [sum(getindex.(possible_machines[inds], [m])) for m in all_m]

	a,b = extrema(A/length(inds))

	state = (b<1.5a) ? 1 : nothing

	return next_m(possible_machines, inds, all_m, min_m, outside_m, state)

end 

function next_m(possible_machines, inds,
								all_m::AbstractVector{<:AbstractString}, 
								min_m::Int, outside_m::Int,
								state
								)::Tuple{<:Any,<:Any}

	next_m(all_m, length(all_m)-outside_m:-1:min_m, state)

end 

function next_m(all_m::AbstractVector{<:AbstractString}, 
								ns::AbstractVector{Int}, i::Int,
								)::Tuple{<:Any,<:Any}


	i>length(ns) ? (nothing,nothing) : (all_m[1:ns[i]], i+1)

end 


function next_m(all_m::AbstractVector{<:AbstractString}, 
								ns::AbstractVector{Int}, state::Union{Nothing,Tuple}
								)::Tuple{<:Any,<:Any}

	it = Iterators.flatten((Combinatorics.combinations(all_m, n) for n=ns))

	S = isnothing(state) ? iterate(it) : iterate(it, state)

	return isnothing(S) ? (nothing,nothing) : S

end 



function Base.iterate(q::IterSplitWS, state=(q.i, nothing))
	
	i, pos, = state 

	i in q.i:length(q.sectors) || return nothing 

	sector = q.sectors[i]

	pos = Assign_Value(pos, i==q.i ? length(sector) : 1)
	
	pos in axes(sector,1) || return iterate(q,(i+1,1))

	stop = sector[pos] 

	stop<q.start && return iterate(q,(i+1,1))


	leftout_ms = i==q.i ? length(sector)-pos : 0


	m,state_comb = next_m(q.possible_machines, q.start:stop,
												q.m, q.min_ms, leftout_ms, state[3:end]...)



	isnothing(m) && return iterate(q, (i, pos + 1 - 2(i==q.i)))
	
	C = cost(q.M, q.possible_machines, q.start:stop, m)


	if length(m)==length(q.m)-leftout_ms && !q.cconstraint(C)

		return i>q.i ? nothing : iterate(q, (i,pos-1))

	end  


	return (i, q.start:stop, m, C), (i, pos, state_comb)

end 



function possible_extensions_(data::Tuple{<:Any,<:Any,<:Any}, 
															candidate::AbstractVector)

	isempty(candidate) && return IterSplitWS(data, 
																	 available_machines(data, candidate, 1))

	prev_i, prev_range, prev_m, costc = candidate[end]

	start = prev_range[end] + 1 

	start > data[2][end][end] && return []

	return IterSplitWS(data,
						 available_machines(data, candidate, start),
						 1,
						 start,
						 prev_i + !in(start, data[2][prev_i]),
						 costc,
						 )
end 


function possible_extensions(data, candidate) 

	(vcat(candidate, [x]) for x in possible_extensions_(data, candidate))

end 



#===========================================================================#
#
# filter candidates/solutions 
#
#---------------------------------------------------------------------------#





function promising_candidate(data, candidate)

	isempty(candidate) && return true 

	#sum(length,getindex.(candidate,3)) > 6 && return false  

#	for (i, sector, machines, cost) in candidate 
#
#		length(machines)>length(sector) && return false 
#
#	end 

	C = getindex.(candidate,4)

	return cost_constraint(C) && cost_constraint2(C)

end  



function accept_sol((possible_machines, sectors, M), candidate)::Bool

	isempty(candidate) && return false 

	sectors[end][end]==candidate[end][2][end] || return false 


	most_powerful = collect(keys(possible_machines[end]))[1:length(candidate)] 

	used_machines = union(getindex.(candidate, 3)...)

	return issubset(most_powerful, used_machines)

end 



function nr_procs_launched(possible_machines, candidate)::Int

	#Vector{Tuple{Int,UnitRange,Vector{String},Float64}}

	n = 0 

	for (i, sector, machines, cost) in candidate 

		for m in machines

			n += minimum(possible_machines[s][m] for s in sector)

		end 

	end 

	return n 

end 


function duplicate_solution(A, B)::Bool

	length(A)==length(B) || return false 

	for ((ia,sa,ma,ca),(ib,sb,mb,cb)) in zip(A,B)

	#	ia==ib || return false 

		length(ma)==length(mb) || return false 

		L = 0.9length(intersect(sa,sb))

		(length(sa)>L && length(sb)>L) || return false 

		(Set(ma)==Set(mb) || isapprox(ca,cb,rtol=0.2)) || return false 

	end 

	return true 

end 


function output((machines, s, M), solutions, candidate)

	any(s->duplicate_solution(s, candidate), solutions) && return true

#	return push!(solutions, candidate)

	length(solutions)<10 && return push!(solutions, candidate)

	n = nr_procs_launched(machines, candidate) 

	for (i,s) in enumerate(solutions) 

		nr_procs_launched(machines, s)<n || continue 

		return setindex!(solutions, candidate, i)

	end 


end 

function available_machines((possible_machines, sectors, M), candidate, start::Int)

	setdiff(collect(keys(possible_machines[start])), getindex.(candidate,3)...)

end 


function split_in_tiers(M::AbstractVector{<:Real},
												possible_machines::AbstractVector{<:AbstractDict},
											 )#::Vector{Tuple{Vector{String},Vector{Int}}}

	s1 = [issorted(M,rev=rev) for rev=[true,false]]
	
	s2 = [issorted(possible_machines, by=length, rev=rev) for rev=[true,false]]

	@assert any(s1) & any(s2)

	#	difficulty: M descending and and nr machines increasing 
	


	split_in_tiers_(s1[2] ? reverse(M) : M,
									s2[1] ? reverse(possible_machines) : possible_machines,
									)


end 


function split_in_tiers_(M::AbstractVector{<:Real},
												 possible_machines::AbstractVector{<:AbstractDict},
												 )::Vector{Vector{Tuple{Vector{Pair{String,Int}},UnitRange{Int}}}}
	
	#Min/max workstations per level and in total
#Tol_cost_gradient
#Min max tiers
#Strictness duplicate
#Nr sols, index sol, sort sola by

	solutions = Backtracking((possible_machines,
														IdentifySectors(keys.(possible_machines)),
														M),
													 possible_extensions,
													 promising_candidate,
													 accept_sol,
													 output,
													 data->[])

	inds = sortperm(solutions, by=s->nr_procs_launched(possible_machines,s))


	return map(reverse(inds)) do j

		map(reverse(solutions[j])) do (i, sector, ws, c)

			([w=>minimum(possible_machines[s][w] for s in sector) for w in ws],
			 
			 UnitRange(extrema(1+length(M).-sector)...)
			 )

		end 

	end 


end 




#############################################################################
end # split workstations 
#############################################################################





#############################################################################

end
