module Taylor
#############################################################################

import Base: == 


import ..LA, ..Utils 


#===========================================================================#
#
# Basic structures 
#
#---------------------------------------------------------------------------#


struct FieldPlaceholder{Rank}

	Name::String 

end 

struct Index{Rank}

	I::NTuple{Rank,Int}

end 

struct Product{Rank}

	Arg::FieldPlaceholder{Rank}

	Weight::Union{Float64,ComplexF64}

	Inds::Vector{Index{Rank}}

end 


struct MixedProduct{NrArgs}

	Weight::Union{Float64,ComplexF64}

	Factors::NTuple{NrArgs,Product}

	function MixedProduct{N}(w::Union{Float64,ComplexF64},
													ps::NTuple{M,Product}
													)::MixedProduct{N} where {N,M}

		@assert N==M && allunique(fieldargs.(ps))

		return new{N}(w, ps)

	end 
 
end   


struct Scalar{NrArgs,NrTerms}

	Args::NTuple{NrArgs,FieldPlaceholder}

	FieldDistrib::NTuple{NrTerms,Vector{Int}} 

	Weight::Union{Float64,ComplexF64} 

	Terms::NTuple{NrTerms,MixedProduct}

end 


struct Tensor{D}
	
	Weight::Union{Float64,ComplexF64}
	
	Dimensions::NTuple{D,Int}

	Inds::Vector{NTuple{D,Index}}

	Components::Vector{Scalar}

end 


#===========================================================================#
#
# the action of the structs when used as functions 
#
#---------------------------------------------------------------------------#


function (i::Index{R})(field::AbstractArray{T,R}
											 )::T where {R,T}#<:Union{Float64,ComplexF64}}

	field[i.I...]

end  

function (f::FieldPlaceholder{R})(field::A
																	)::A where {R,
																						 T<:Union{Float64,ComplexF64},
																						 A<:AbstractArray{T,R}}
	field

end 

function (p::Product)(field::AbstractArray{T}
												)::Union{Float64,ComplexF64} where T

	small_weight(p) && return 0.0

	out::promote_type(typeof(p.Weight),T) = p.Weight

#	fieldargs(p)(field)  # check shape 

	for I in each_fieldfactor(p) 
		
		out *= I(field)

	end 


	return out 

#	return prod(I(field) for I in each_fieldfactor(p);init=p.Weight)
#	return prod(field)

end 

function (P::MixedProduct{N})(fields::Vararg{<:AbstractArray,N}
																 )::Union{Float64,ComplexF64} where N 

	small_weight(P) && return 0.0

	out = P.Weight 

	for (p,f) in zip(parts(P),fields) 

		out *= p(f)

	end 

	return out 

end 




function (S::Scalar{N})(fields::Vararg{<:AbstractArray, N}
													 )::Union{Float64,ComplexF64} where N 

	small_weight(S.Weight) && return 0.0 

	s = 0.0 

	for (p,f) in zip(parts(S), fieldargs_distrib(S))

		s += p(fields[f]...)

	end 

	return s*S.Weight 

end 

function (PPS::Union{Product,MixedProduct,Scalar})(fields::Tuple{Vararg{<:AbstractArray}})::Number 

	PPS(fields...)

end




#===========================================================================#
#
#	zeros 
# 
#---------------------------------------------------------------------------#


function Base.zero(p::Product{R})::Product{R} where R

	Product{R}(fieldargs(p), 0.0, Vector{Index{R}}(undef,0))

end 

function Base.zero(P::T)::T where T<:MixedProduct

	T(0.0, zero.(parts(P)))

end 


function zero_scalar(args::NTuple{Na,FieldPlaceholder},
										 fd::NTuple{Nt, AbstractVector{Int}},
										 terms::NTuple{Nt, MixedProduct}
										 )::Scalar{Na,1} where {Na,Nt}

	i = argmin(length.(terms))

	return Scalar{Na,1}(args, fd[[i]], 0.0, tuple(zero(terms[i])))

end 


function Base.zero(S::Scalar{Na})::Scalar{Na,1} where Na

	zero_scalar(fieldargs(S), fieldargs_distrib(S), parts(S))

end 






#===========================================================================#
#
# field rank 
#
#---------------------------------------------------------------------------#

function field_rank(I::NTuple{N,Int})::Int where N

	N 

end 

function field_rank(I::Index{R})::Int where R

	R

end  

function field_rank(I::AbstractVector{<:Union{Index,Tuple{Vararg{Int}}}}
									 )::Int 

	r = field_rank(I[1])

	for i in 2:length(I)

		@assert r == field_rank(I[i])

	end 

	return r 

end 

nr_fields(I::AbstractVector{<:Index})::Int = length(I)



function each_fieldfactor(I::T)::T where T<:AbstractVector{<:Index}

	I

end 

function field_rank(f::FieldPlaceholder{R})::Int where R 

	R

end 

function field_rank(p::Product{R})::Int where R
	
	R

end  


function field_rank(c::MixedProduct{N})::NTuple{N,Int} where N
	
	map(field_rank, parts(c))

end 

function nr_fields(c::MixedProduct{N})::NTuple{N,Int} where N
	
	map(nr_fields, parts(c))

end 


nr_fields(c::MixedProduct, field::Int)::Int = nr_fields(parts(c,field))

field_rank(c::MixedProduct, field::Int)::Int = field_rank(parts(c,field))



enumerate_fieldfactors = enumerate∘each_fieldfactor 
nr_fields(p::Product)::Int = nr_fields(p.Inds)  



#===========================================================================#
# 
# fieldargs 
#
#---------------------------------------------------------------------------#



function fieldargs(p::Product{R})::FieldPlaceholder{R} where R 

	p.Arg 

end 


function fieldargs(P::MixedProduct{N})::NTuple{N,FieldPlaceholder} where N
	
	fieldargs.(parts(P))

end 

function fieldargs(S::Scalar{Na})::NTuple{Na,FieldPlaceholder} where Na

	S.Args 
	
end 

function fieldargs(S::Scalar{Na}, i::Int)::FieldPlaceholder where Na

	@assert 1<=i<=Na 

	return fieldargs(S)[i]

end 

#===========================================================================#
#
# Defining what the "parts" of each object mean 
#
#---------------------------------------------------------------------------#


function parts(p::Product{R})::AbstractVector{Index{R}} where R

	each_fieldfactor(p)

end 


function parts(P::MixedProduct{N})::NTuple{N,Product} where N
	
	P.Factors 

end 

function parts(P::MixedProduct{N}, i::Int)::Product where N 
	
	@assert 1<=i<=N 
	
	return P.Factors[i] 

end 

function parts(S::Scalar{Na,Nt})::NTuple{Nt,MixedProduct} where {Na,Nt}
	
	S.Terms 
 
end  

function parts(S::Scalar{Na,Nt}, i::Int)::MixedProduct where {Na,Nt}

	@assert 1<=i<=Nt 

	return S.Terms[i]

end 

function parts(T::Tensor)::Vector{Scalar}

	T.Components

end 

function parts(T::Tensor, i::Int)::Scalar 

	@assert 1<=i<=length(T.Inds) 

	return T.Components[i]

end 

function parts(T::Tensor, ::Nothing)::Scalar 

	zero(parts(T,1))

end 




Base.length(a::Union{Product,MixedProduct,Scalar,Tensor})::Int = length(parts(a))

enum_parts = enumerate∘parts  

#===========================================================================#
#
# field distribution 
#
#---------------------------------------------------------------------------#


function fieldargs_distrib(S::Scalar{Na,Nt}
													 )::NTuple{Nt,Vector{Int}} where {Na,Nt}

	S.FieldDistrib 

end 

function fieldargs_distrib(S::Scalar{Na,Nt}, i::Int
													)::Vector{Int} where {Na,Nt}

	@assert 1<=i<=Nt

	return fieldargs_distrib(S)[i]

end  


function fieldargs_distrib(S::Scalar{Na,Nt}, I::AbstractVector{Int}
													 )::Tuple{Vararg{Vector{Int}}} where {Na,Nt}

	for i in I 
		
		@assert 1<=i<=Nt

	end 

	return fieldargs_distrib(S)[I]

end 

function fieldargs_distrib(fph::Union{AbstractVector{<:FieldPlaceholder},
																			Tuple{Vararg{<:FieldPlaceholder}}}, 
													 args...)

	fieldargs_distrib([argnames(a) for a in fph], args...)

end


function fieldargs_distrib(unique_names::AbstractVector{<:AbstractString},
													 item::AbstractString,
													 )::Vector{Int} 

	fieldargs_distrib(unique_names, [item])

end 

function fieldargs_distrib(unique_names::AbstractVector{<:AbstractString},
													 item::Union{AbstractVector{<:AbstractString},
																			 Tuple{Vararg{<:AbstractString}}},
													 )::Vector{Int} 
	
	i = indexin(item, unique_names)

	@assert all(!isnothing, i) "Some fields not found"

	return Vector{Int}(i) 

end  



function fieldargs_distrib(unique_names::AbstractVector{<:AbstractString},
													 v::Tuple{Vararg{AbstractVector{<:AbstractString}}},
													 )::Tuple{Vararg{Vector{Int}}}

	Tuple(fieldargs_distrib(unique_names, item) for item in v)

end  


function fieldargs_distrib(unique_names::AbstractVector{<:AbstractString},
													 p::Union{Product,MixedProduct},
													 )::Vector{Int}

	fieldargs_distrib(unique_names, argnames(p))

end 

function fieldargs_distrib(unique_names::AbstractVector{<:AbstractString},
													 v::Union{Tuple{Vararg{<:Union{<:Product,MixedProduct}}},
																		AbstractVector{<:Union{<:Product,MixedProduct}}},
													 )::Tuple{Vararg{Vector{Int}}}

	Tuple(fieldargs_distrib(unique_names, item) for item in v)

end  


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#




#===========================================================================#
#
# names for the fields 
#
#---------------------------------------------------------------------------#



argnames(s::AbstractString)::String = s 

argnames(c::Char)::String = string(c) 

argnames(rank::Int)::String = string('A'+rank-1) 

function argnames(M::AbstractVector{<:Index})::String# where N 

	argnames(field_rank(M))

end 


argnames(s::Symbol)::String = string(s)

argnames(f::FieldPlaceholder)::String = f.Name



argnames(P::Product)::Tuple{String} = tuple(argnames(fieldargs(P)))

function argnames(P::T)::NTuple{N,String} where {N,T<:Union{MixedProduct{N},Scalar{N}}}

	argnames.(fieldargs(P))

end 


function argnames(P::Union{MixedProduct,Scalar}, i::Int)::Tuple{Vararg{String}}

	argnames(parts(P,i))

end 




#===========================================================================#
#
# Utils for Index struct 
#
#---------------------------------------------------------------------------#


function disregard_fieldfactors(I::AbstractVector{<:Index},
																i::Union{Int,AbstractVector{Int}},
																)::AbstractVector{<:Index}

	select_fieldfactors(I, setdiff(1:length(I), i))

end 


function select_fieldfactor(I::AbstractVector{<:Index},
														 i::Int,
											 )::Index 

	I[i]

end 

function select_fieldfactors(I::AbstractVector{<:Index},
														 i::Union{Int,AbstractVector{Int}},
														 )::AbstractVector{<:Index}

	view(I, vcat(i))

end 



function split_fieldfactor(I::AbstractVector{<:Index},
														i::Int
														)::Tuple{Index,AbstractVector{<:Index}}

	(select_fieldfactor(I,i), disregard_fieldfactors(I,i))

end 



function each_fieldfactor(p::Product{R})::Vector{Index{R}} where R 
	
	each_fieldfactor(p.Inds)

end 


function disregard_fieldfactors(p::Product{R}, args...
																)::AbstractVector{Index{R}} where R 
	
	disregard_fieldfactors(p.Inds, args...)

end 


function select_fieldfactors(p::Product{R}, args...
														 )::AbstractVector{Index{R}} where R

	select_fieldfactors(p.Inds, args...)

end 

function select_fieldfactor(p::Product{N}, args...)::Index{N} where N
	
	select_fieldfactor(p.Inds, args...)

end 

function split_fieldfactor(P::Product{R}, i::Int
													 )::Tuple{Index{R}, AbstractVector{Index{R}}
																		} where R

	split_fieldfactor(P.Inds, i)

end 



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function Base.in(i::Index{N}, p::Product{R})::Bool where {N,R}

	N==R && in(i,each_fieldfactor(p))

end  

function Base.in(f::FieldPlaceholder, X::Union{MixedProduct,Scalar}
								)::Bool 

	f in fieldargs(X)

end 

function Base.findfirst(i::Index{N}, p::Product{R}
												)::Union{Int,Nothing} where {N,R}

	if N==R 

		for (i0,I0) in enumerate_fieldfactors(p)
	
			i==I0 && return i0 
	
		end 

	end 

	return nothing 

end 

function Base.findall(
											i::Index{N}, p::Product{R}
											)::Vector{Int} where {N,R}

	N!=R ? Int[] : [i0 for (i0,I0) in enumerate_fieldfactors(p) if i==I0] 

end  

function Base.findall(f::FieldPlaceholder, X::Union{MixedProduct,Scalar}
											)::Vector{Int}

	findall(==(f), fieldargs(X))

end 

function Base.findfirst(f::FieldPlaceholder, 
												X::Union{MixedProduct,Scalar}
												)::Union{Nothing,Int}

	findfirst(==(f), fieldargs(X))

end 

function Base.findfirst(fn::AbstractString,
												X::Union{MixedProduct,Scalar}
												)::Union{Nothing,Int}

	findfirst(==(fn), argnames(X))

end 

function Base.findfirst(I::NTuple{N,Index}, T::Tensor{D}
												)::Union{Nothing,Int} where {N,D}

	N==D ? findfirst(==(I), T.Inds) : nothing 

end 

function Base.join(I::NTuple{N,Index}, i::Index)::NTuple{N+1,Index} where N 

	(I..., i)

end  

#===========================================================================#
#
# equality and proportionality 
#
#---------------------------------------------------------------------------#

function proportional(p::P, q::Q
											 )::Bool where {T<:Union{Product, 
																							MixedProduct, 
																							Scalar},P<:T,Q<:T}

	same_fields(p,q) && Utils.has_disjoint_pairs(proportional, parts(p), parts(q))

end 


function proportional(p::Union{Product,MixedProduct,Scalar})

	function prop(q::Union{Product,MixedProduct,Scalar})::Bool
		
		proportional(p,q)

	end 

end  

function same_inds(p::Product, q::Product)::Bool 

	Utils.has_disjoint_pairs(==, parts(p), parts(q))

end 

function ==(i1::Index, i2::Index)::Bool# where T<:Tuple{Vararg{Int}}

	field_rank(i1)==field_rank(i2) && i1.I==i2.I

end  



proportional(i1::Index, i2::Index)::Bool = i1==i2 


function equal_prods(t1::T, t2::T, 
										 W1::Union{Real,Complex},
										 W2::Union{Real,Complex}
										 )::Bool where T<:Union{Product,MixedProduct}

	proportional(t1,t2) && same_weight(weight(t1)*W1, weight(t2)*W2)

end  


function equal_prods(S1::Scalar{Na1,Nt}, S2::Scalar{Na2,Nt}
										 )::Bool where {Na1,Na2,Nt}
	
	Utils.has_disjoint_pairs(equal_prods, parts(S1), parts(S2), S1.Weight, S2.Weight)

end 

function equal_prods(S1::Scalar, S2::Scalar)::Bool 

	@assert  length(S1)!=length(S2)

	return false  

end 



function ==(P::T1, Q::T2
					 )::Bool where {T<:Union{Product,MixedProduct},T1<:T,T2<:T}

	if same_fields(P,Q) && same_weight(P,Q)
		
		return small_weight(P) ? true : proportional(P, Q) 
		
	end 
	
	return false 

end 



function ==(S1::Scalar, S2::Scalar)::Bool # assuming unique terms 

	same_fields(S1,S2) || return false 

	small_weight(S1.Weight) && small_weight(S2.Weight) && return true 

	return equal_prods(S1, S2)  

end 



#===========================================================================#
#
# unique 
#
#---------------------------------------------------------------------------#

function count_unique(p::Product)::Vector{NTuple{2,Int}}#Product}

	U,I = Utils.find_degen(select_fieldfactor, nr_fields(p), p)

	return [(u,length(i)) for (u,i) in zip(U,I)]

end 


#===========================================================================#
#
# 	comparing fields 
#
#---------------------------------------------------------------------------#



function ==(f1::FieldPlaceholder{N1}, 
						f2::FieldPlaceholder{N2})::Bool where {N1,N2}
	
	if argnames(f1)==argnames(f2)

		@assert N1==N2 "Same name used for different fields"

		return true 

	end 

	return false 

end  





function same_fields(p::P, q::Q)::Bool where {T<:Union{Product,MixedProduct,Scalar},P<:T,Q<:T}

	fp::Union{Tuple,FieldPlaceholder} = fieldargs(p)

	fq::Union{Tuple,FieldPlaceholder} = fieldargs(q) 

	return Utils.has_disjoint_pairs(==,
														fp isa Tuple ? fp : tuple(fp),
														fq isa Tuple ? fq : tuple(fq)
														)

end 


function same_fields(p::Union{Product,MixedProduct,Scalar})::Function

	function same_fields_(q::Union{Product,MixedProduct,Scalar})::Bool

		same_fields(p,q)

	end 

end 





#===========================================================================#
#
# weight 
#
#---------------------------------------------------------------------------#

small_weight(w::Number)::Bool = abs2(w)<1e-20

same_weight(w1::Number, w2::Number)::Bool = small_weight(w1-w2)  


weight(p::Product)::Union{Float64,ComplexF64} = p.Weight

weight(P::MixedProduct)::Union{Float64,ComplexF64} = P.Weight*prod(weight, parts(P))

weight(S::Scalar, i::Int)::Union{Float64,ComplexF64} = weight(parts(S,i))


function small_weight(p::Union{Product,MixedProduct})::Bool 
	
	small_weight(weight(p))

end 


function same_weight(p::Union{Product,MixedProduct},
										 q::Union{Product,MixedProduct})::Bool 

	same_weight(weight(p), weight(q))

end 





#===========================================================================#
#
# "Index" constructors 
#
#---------------------------------------------------------------------------#

function parse_inds(inds...
										)::Vector{<:Index}

	rank = map(inds) do i
		
		@assert i isa Union{Int, Tuple{Vararg{Int}}, AbstractVector{Int}}   
	
		return length(i)

	end 

	@assert length(unique(rank))==1 

	R = rank[1] 

	return [Index{R}(Tuple(getindex(i,r) for r=1:R)) for i in inds] 

end  

Index(I::Index)::Index = I 
Index(i::Int)::Index{1} = Index(tuple(i))
Index(v::AbstractVector{Int})::Index = Index(Tuple(v))

Index(I::CartesianIndex)::Index = Index(Tuple(I))



#===========================================================================#
#
# "Product" constructors 
#
#---------------------------------------------------------------------------#



function Product(f::FieldPlaceholder{N},
										w::Union{Float64,ComplexF64},
										M::AbstractVector{<:Union{Tuple,<:Index}},#NTuple{N,Int}},
										)::Product where N

	Inds = Vector{Index{N}}(undef,length(M))

	for (i,m) in enumerate(M)

		Inds[i] = Index(m)

	end 

	return Product{N}(f, w, Inds)

end 

function Product(name::Union{AbstractString,Symbol,Char},
										w::Union{Float64,ComplexF64},
										M::AbstractVector{<:Union{Tuple,<:Index}},#NTuple{N,Int}},
										)::Product #where N

	Product(FieldPlaceholder{field_rank(M)}(argnames(name)), w, M)

end 


function Product(w::Union{Float64,ComplexF64},
										M::AbstractVector{<:Union{Tuple,<:Index}},#NTuple{N,Int}},
										)::Product #where N

	Product(argnames(M), w, M)

end 


function Product(name::Union{Char,AbstractString,Symbol},
										M::AbstractVector{<:Union{Tuple,<:Index}}
										)::Product 

	Product(argnames(name), 1.0, M)

end 


function Product(M::AbstractVector{<:Union{Tuple,<:Index}}
								)::Product 

	Product(1.0, M)

end 



function Product(name::Union{Char,AbstractString,Symbol},
										w::Union{Float64,ComplexF64},
										ind1::Union{Int,Tuple{Vararg{Int}},AbstractVector{Int}},
										inds...)::Product 

	Product(argnames(name), w, parse_inds(ind1, inds...))

end 

function Product(w::Union{Float64,ComplexF64}, 
										ind1::Union{Int,Tuple{Vararg{Int}},AbstractVector{Int}},
										inds...
										)::Product 

	Product(w, parse_inds(ind1, inds...))

end 



function Product(ind1::Union{Int,Tuple{Vararg{Int}},AbstractVector{Int}},
										inds...)::Product 

	Product(parse_inds(ind1, inds...))

end  




function Product(name::Union{Char,AbstractString,Symbol},
										ind1::Union{Int,Tuple{Vararg{Int}},AbstractVector{Int}}, 
										inds...)::Product 

	Product(name, parse_inds(ind1, inds...))

end 



function Product_sameField(ps::Vararg{Product})::Product 

	Product_sameField(ps)

end 
								
function Product_sameField(ps::Union{Tuple{Vararg{<:Product}},
																				AbstractVector{<:Product}}
															)::Product 

	length(ps)==1 && return first(ps)

	f = fieldargs(first(ps))

	@assert all(==(f)∘fieldargs, ps) "All prods must take the same arg"

	return Product(f, 
										prod(weight, ps; init=1.0),
										mapreduce(parts,vcat,ps))


end 





function parse_inds_2(w::T, inds_or_prods...
											)::Tuple{T,Vector{Product}
															 } where T<:Union{Float64,ComplexF64}

	(w, map(1:length(inds_or_prods)) do i 

		 q = inds_or_prods[i]

		 q isa Product && return q 

		 Utils.isListLen(q, Product, 1) && return only(q)

		 return Product(q...)

	end)

end 

function parse_inds_2(arg1::Union{Product, Int, Index,
																	Tuple{Vararg{Int}}, AbstractVector{Int}},
											args...
										 )::Tuple{Float64,Vector{Product}}

	parse_inds_2(1.0, arg1, args...)

end  




#===========================================================================#
#
# constructors for Mixed Products 
#
#---------------------------------------------------------------------------#


function MixedProduct_(w::Union{Float64,ComplexF64},
													ps::NTuple{N,Product}
												 )::MixedProduct{N} where N 

	@assert allunique(fieldargs.(ps)) "The args cannot be used multiple times"

	return MixedProduct{N}(w, ps)

end  



function MixedProduct(fns::AbstractVector{<:Union{Char,AbstractString,Symbol}},
												 w::Union{Float64,ComplexF64},
												 ps::Union{AbstractVector{<:Product},
																	 Tuple{Vararg{<:Product}}}
												 )::MixedProduct

	unique_field_names = unique(argnames.(fns)) 

	I = only.(fieldargs_distrib(unique_field_names, ps)) 

#	I[2]= 5 means the second product uses the fifth unique field 


	items, positions = Utils.Unique(I, sorted=true, inds=:all) 

	@assert items==1:length(unique_field_names)


	return MixedProduct_(w, Tuple(Product_sameField(ps[i]) for i in positions)
												 )
end 




function MixedProduct(w::Union{Float64,ComplexF64}, 
												 v::Union{AbstractVector{<:Product},
																	Tuple{Vararg{<:Product}}}
												)::MixedProduct

	MixedProduct([only(argnames(p)) for p in v], w, collect(v)) 

end 

function MixedProduct(arg1::Union{Float64, ComplexF64, Product}, 
												 args...
												)::MixedProduct  

	MixedProduct(parse_inds_2(arg1, args...)...)

end 

function MixedProduct(args::Tuple{Vararg{Product}})::MixedProduct

	MixedProduct(parse_inds_2(args...)...)

end 

function MixedProduct(fns::Tuple,
												 arg1::Union{Float64,ComplexF64,Product},
												 args...)::MixedProduct

	MixedProduct([argnames(fn) for fn in fns], 
									parse_inds_2(arg1, args...)...)

end 

function MixedProduct(fns::AbstractVector{<:Union{Char,AbstractString,Symbol}},
												 arg1::T1,
												 arg2::T2,
												 args...)::MixedProduct where {T<:Union{Float64,ComplexF64, Product, Int, Tuple{Vararg{Int}}, Index, AbstractVector{Int}},T1<:T,T2<:T}

	MixedProduct(fns, parse_inds_2(arg1, arg2, args...)...)
	
end 





#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#




function MixedProduct_sameFields(P::MixedProduct, 
																		p1::Product, 
																		i1::Int,
																		args...
												 )::MixedProduct

	MixedProduct_sameFields_(P, p1, i1, args...)

end 

function MixedProduct_sameFields_(P::MixedProduct, 
																		 args::Vararg{T,N} where T
																		 )::MixedProduct where N 

	@assert iseven(N)

	inds, locs = Utils.Unique([args[2i] for i=1:div(N,2)]; inds=:all)

	return MixedProduct_(P.Weight, ntuple(length(P)) do i 

		k = findfirst(isequal(i),inds)
		
		isnothing(k) && return parts(P,i)

		return Product_sameField(parts(P,i), (args[2j-1] for j in locs[k])...)

	end)


end 



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function replace_component(P::T,
													 new_factor::Product,
													 location::Int)::T where T<:MixedProduct  

	MixedProduct_(P.Weight, ntuple(length(P)) do i 

										 i==location ? new_factor : parts(P,i)

										end)

end 





#===========================================================================#
#
# Scalar constructors 
#
#---------------------------------------------------------------------------#



function Scalar_(args::NTuple{Na,FieldPlaceholder},

									 fd::NTuple{Nt, AbstractVector{Int}},

									 w::Number,

									 terms::NTuple{Nt, MixedProduct},

									 )::Scalar where {Na,Nt}

	categories = cumulate_categories(terms) 
	
	u_inds = unique_cumul_categ(categories)


	if isempty(u_inds) || small_weight(w)

		return zero_scalar(args, fd, terms)

	elseif categories==1:length(terms) 

		return Scalar{Na,Nt}(args, fd, w, terms)
		
	else 

		return Scalar{Na,length(u_inds)}(args,
																	Tuple(fd[i[1]] for i in u_inds),
																	w,
																	cumulate_(terms, u_inds)
																	)
	end 

end 

function Scalar(w::Number,
									 unique_fields::Union{AbstractVector{<:FieldPlaceholder},
																			 Tuple{Vararg{<:FieldPlaceholder}}},
									 terms::Union{AbstractVector{<:MixedProduct},
																Tuple{Vararg{<:MixedProduct}}}
																)::Scalar 
	

	field_distrib = fieldargs_distrib(unique_fields, terms)


	for index_unique_field in 1:length(unique_fields)
		
#		map(field_distrib) do fields_of_term_j 
		for fields_of_term_j in field_distrib

			pos_field_i_in_term_j = findall(fields_of_term_j.==index_unique_field) 

			@assert length(pos_field_i_in_term_j)<2 "No field can be used twice"

		end 

	end  

	return Scalar_(Tuple(unique_fields), field_distrib, w, Tuple(terms))

end 



function Scalar(w::Number,
									 terms::Union{AbstractVector{<:MixedProduct},
																Tuple{Vararg{<:MixedProduct}}},
																)::Scalar   

	@assert !isempty(terms) #&& return zero(terms)

	return Scalar(w, unique(mapreduce(fieldargs, Utils.tuplejoin, terms)), terms)

end 


function Scalar(arg1::Union{Tuple,AbstractVector,MixedProduct}, 
									 args...)::Scalar 
	
	Scalar(1.0, arg1, args...)

end 


function Scalar(w::Number, args::Vararg{MixedProduct})::Scalar

	Scalar(w, args)	#	 [p for p in args])

end 



#===========================================================================#
#
# Tensor constructors 
#
#---------------------------------------------------------------------------#


function Tensor(w::Number, s::NTuple{D,Int}, S::Scalar
									)::Tensor{D} where D 

	Tensor{D}(w, s, [Tuple(Index(fill(1,r)) for r in s)], [S])

end 

function Tensor(S::Scalar)::Tensor{0}

	Tensor(1.0, (), S)

end 

#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


														 
function Tensor_fromDeriv(S::Scalar, 
														 fields::Vararg{Any,D}
														)::Tensor{D} where D

	Tensor_fromDeriv(Tensor(S), fields...)

end 

function Tensor_fromDeriv(T::Tensor{D},
														 field::Union{AbstractString,FieldPlaceholder},
														 fields::Vararg{Any,N}
														 )::Tensor{D+N+1} where {D,N}

	Tensor_fromDeriv(Tensor_fromDeriv(T,field), fields...)

end 


function Tensor_fromDeriv(T::Tensor{D}, field::AbstractString
														 )::Tensor{D+1} where D

	for S in parts(T) 

		i = findfirst(field, S) 

		i isa Int && return Tensor_fromDeriv(T,  fieldargs(S)[i])

		for P in parts(S) 

			j = findfirst(field, P) 

			j isa Int && return Tensor_fromDeriv(T, fieldargs(P)[j])

		end 

	end 

	error("The rank of the field not known")

end 



function Tensor_fromDeriv(T::Tensor{D}, field::FieldPlaceholder{R}
														 )::Tensor{D+1} where {D,R}

	inds, scalars = derivatives(T, field)
	
	s = (size(T)..., R)

	isempty(inds) && return Tensor(0.0, s, parts(T,nothing))

	return Tensor(1.0, s, inds, scalars)

end 





#===========================================================================#
#
# tensor properties 
#
#---------------------------------------------------------------------------#





function Base.ndims(T::Tensor{D})::Int where D 
	
	D

end 

function Base.size(T::Tensor{D})::NTuple{D,Int} where D 

	T.Dimensions

end 

function Base.size(T::Tensor{D}, i::Int)::Int where D

	@assert 1<=i<=D

	size(T)[i]

end 

#Base.length(T::Tensor)::Int = length(T.Inds)

function Base.getindex(T::Tensor, I::Vararg{Index})::Scalar 

	T[I]

end 

function Base.getindex(T::Tensor, I::Tuple{Vararg{Index}})::Scalar 

	parts(T,findfirst(I,T))

end 


function Base.iterate(T::Tensor, state::Int=1)

	state > length(T) && return nothing 

	return (T.Inds[state],parts(T,state)), state+1

end 




#===========================================================================#
#
# derivatives 
#
#---------------------------------------------------------------------------#

function derivative(S::Scalar{N}, field_index::Int, args...
										)::Scalar{N} where N 

	derivative(S, fieldargs(S, field_index), args...)

end 


function derivative(S::Scalar{N}, name::AbstractString, args...
										)::Scalar{N} where N 

	derivative(S, findfirst(name, S), args...)

end 


function derivative(S::Scalar{N}, 
										field::FieldPlaceholder{R},
										comp::Index{R}
										)::Scalar{N} where {R,N}

	field in S || return zero(S)

	return Scalar_(fieldargs(S), fieldargs_distrib(S), S.Weight,
										Tuple(derivative(P, field, comp) for P in parts(S)))

end 

function flatderiv_get_I(d::AbstractVector, K::AbstractVector, k::Int)::Index 

	i = flatderiv_get_i(K, k)
	j = flatderiv_get_j(K, k)

	return d[i][1][j]

end 

function flatderiv_get_X(d::AbstractVector, K::AbstractVector, k::Int
								)

	i = flatderiv_get_i(K, k)
	j = flatderiv_get_j(K, k)

	return d[i][2][j]

end  


function flatderiv_get_i(K::AbstractVector, k::Int)::Int 

	(i,j) = K[k]

	return i

end 

function flatderiv_get_j(K::AbstractVector, k::Int)::Int 

	(i,j) = K[k]

	return j

end 




function derivatives(S::Scalar{N},
										 field::Union{AbstractString,FieldPlaceholder},
										 weight_factor::Number=1
										 )::Tuple{Vector{<:Index}, Vector{Scalar{N}}} where N

	data = [derivatives(P,field,S.Weight*weight_factor) for P in parts(S)]

	J = [(a,b) for (a,(nz_inds,)) in enumerate(data) for b in 1:length(nz_inds)]


	isempty(J) && return (Vector{Index}(undef, 0), 
												Vector{Scalar{N}}(undef, 0))




	K_i, K_p = Utils.find_degen(flatderiv_get_I, length(J), data, J)

	return ([flatderiv_get_I(data, J, k) for k in K_i],
					[zero(S) + sum(flatderiv_get_X(data, J, k) for k in K) for K in K_p]
					)

end 




function derivative(P::T, y::Int, x::Union{Int,Index}
										)::T where T<:MixedProduct 

	replace_component(P, derivative(parts(P,y), x), y)

end  

function derivative(P::T, ::Nothing, args...)::T where T<:MixedProduct 

	zero(P)

end 

function derivative(P::T, f::FieldPlaceholder, x::Union{Int,Index},
										)::T where T<:MixedProduct

					# unique entries 

	derivative(P, findfirst(f, P), x) 

end  



function derivatives(P::MixedProduct, y::Int,
										 args...
										 )::Tuple{Vector{<:Index}, Vector{MixedProduct}}

	nz_inds, derivs = derivatives(parts(P, y), args...)

	return nz_inds, [replace_component(P, d, y) for d in derivs]

end  

function derivatives(P::MixedProduct, 
										 field::Union{AbstractString,FieldPlaceholder},
										 args...
										 )::Tuple{Vector{<:Index}, Vector{MixedProduct}}

	y = findfirst(field, P)

	if isnothing(y) 

		return (Vector{Index}(undef,0),Vector{MixedProduct}(undef,0))

	else 

		return derivatives(P, y, args...)

	end 

end 




function derivatives(T::Tensor{D}, 
										 field::Union{AbstractString,FieldPlaceholder},
										 )::Tuple{Vector{NTuple{D+1,Index}},
															Vector{Scalar}} where D

	data = [derivatives(S,field,T.Weight) for S in parts(T)]

	K = [(a,b) for (a,(nz_inds,)) in enumerate(data) for b in 1:length(nz_inds)] 
				 
	return ([join(T.Inds[i], data[i][1][j]) for (i,j) in K],
					[data[i][2][j] for (i,j) in K]
					)

end 








function derivative(p::Product, i::Int, 
										weight_factor::Union{Real,Complex}=1
										)::Product 

	@assert 1<=i<=nr_fields(p)  # degeneracies *not* taken into account

	return Product(fieldargs(p),
										p.Weight*weight_factor, 
										disregard_fieldfactors(p,i)) 

end 


function derivative(p::P, I0::Index)::P where P<:Product
	
	occurences = findall(I0, p)
	
	isempty(occurences) && return zero(p) # degeneracies taken into account 
	
	return derivative(p, occurences[1], length(occurences))
	
end  


function derivatives(p::Product{R},
										 weight_factor::Number=1,
										 )::Tuple{Vector{Index{R}}, Vector{Product{R}}
															} where R

	cu = count_unique(p)

	I = Vector{Index{R}}(undef, length(cu))

	P = Vector{Product{R}}(undef, length(cu))


	for (order,(i,nr)) in enumerate(count_unique(p))

		(i_1, i_rest) = split_fieldfactor(p,i)

		I[order] = i_1

		P[order] = Product(fieldargs(p), weight_factor*p.Weight*nr, i_rest)

	end 

	return I,P

end 





#===========================================================================#
#
# cumulate 
#
#---------------------------------------------------------------------------#





function cumulate_(ps::Union{
														 Tuple{Vararg{Union{Product,MixedProduct}}},
														 AbstractVector{<:Union{Product,MixedProduct}}},
									 degen_inds::AbstractVector{Int}
									 )::Union{Product,MixedProduct}

	p = ps[first(degen_inds)]

	length(degen_inds)==1 && return p 
	
	w = W = weight(p)

	for i=2:length(degen_inds)

		W += weight(ps[degen_inds[i]])

	end 

	return p * (W/w)

end 








function cumulate(ps::AbstractVector{T},
									start::Int=1
									)::AbstractVector{T} where T<:Union{Product, MixedProduct}

	for i=start:length(ps)

		js = i+1:length(ps) 

		sim_i = i .+ findall(proportional(ps[i]), view(ps, js)) 

		js_ = setdiff(js, sim_i) 


		wi = weight(ps[i]) 

		W = sum(weight, view(ps, sim_i); init=wi)

		
		if small_weight(wi) 
			
			return cumulate(view(ps, vcat(1:i-1, small_weight(W) ? js_ : js)), i)
		
		else 

			return cumulate(vcat(view(ps, 1:i-1), ps[i]*(W/wi), view(ps,js_)), i+1)

		end 

	end 

	return ps

end 



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function cumulate_categories(ps::Union{AbstractVector{T},Tuple{Vararg{T}}},
														)::Vector{Int} where T<:Union{<:Product,<:MixedProduct}


	rels::AbstractVector{Int}=zeros(Int,length(ps)) # output 


	W = [weight(p) for p in ps]

	seen = small_weight.(W)  # the insignificant components are not considered 


	for (i,p) in enumerate(ps)

		seen[i] && continue # skip if "i" already looked at 

		seen[i] = true  		# "i" is being looked at now

		@assert rels[i]==0  # rels[i] changes only when seen. Error otherwise 


		rels[i] = i					# has significant weight => its own category 


		sim_i = [j for j in i+1:length(ps) if !seen[j]&&proportional(p,ps[j])] 

		isempty(sim_i) && continue 

		

		@assert all(isequal(0), view(rels, sim_i)) # since it hasn't been seen 

		seen[sim_i] .= true  



		if small_weight(W[i]+sum(view(W, sim_i)))
						# weight "i" significant, but the total amounts to zero 
					
			rels[i] = 0 

		else 

			rels[sim_i] .= i 

		end  # rels[...] remain zero if the total weight is small


	end  


	return rels 

end  


function unique_cumul_categ(categories::AbstractVector{Int}
														)::Vector{Vector{Int}}

	u_categ,u_inds = Utils.Unique(categories,inds=:all)

	I = u_categ.>0 

	for (c, (i,)) in zip(view(u_categ,I), view(u_inds,I))

		@assert c==i

	end  

	return view(u_inds, I) 

end 








function cumulate_(ps::AbstractVector{<:T},
									 u_inds::AbstractVector{<:AbstractVector{Int}}
									 )::Vector{T} where T<:Union{Product,MixedProduct}

	[cumulate_(ps, degen_inds) for degen_inds in u_inds]

end 

function cumulate_(ps::Tuple{Vararg{Union{Product,MixedProduct}}},
									 u_inds::AbstractVector{<:AbstractVector{Int}}
									 )::Tuple#{Vararg{T}} where 

	Tuple(cumulate_(ps, degen_inds) for degen_inds in u_inds)

end 





function cumulate_(ps::AbstractVector{T},
									)::AbstractVector{T} where T<:Union{Product,MixedProduct}

	categories = cumulate_categories(ps)

	categories==1:length(ps) && return ps 

	return cumulate_(ps, unique_cumul_categ(categories))
					 
end 

#===========================================================================#
#
# Arithmetics 
#
#---------------------------------------------------------------------------#


function Base.:*(ps::Vararg{Product})::MixedProduct 

	MixedProduct(ps...)

end 

function Base.:*(P::MixedProduct, p::Product)::MixedProduct

	i = findfirst(==(fieldargs(p)), fieldargs(P))

	if isnothing(i) 
		
		@assert all(!=(only(argnames(p))), argnames(P))

		return MixedProduct_(P.Weight, (parts(P)..., p)) 

	else 

		return replace_component(P, Product_sameField(P.Factors[i], p), i)

	end 

end   



function Base.:*(p::Product, P::MixedProduct)::MixedProduct

	i = findfirst(fieldargs(p), P)


	if isnothing(i) 
		
		@assert all(!=(only(argnames(p))), argnames(P))

		return MixedProduct_(P.Weight, (p, parts(P)...)) 

	else 

		return replace_component(P, Product_sameField(p, parts(P, i)), i)

	end 

end  



function Base.:*(w::Number, p::T)::T where T<:Union{Product,
																										MixedProduct,
																										Scalar}
	same_weight(1,w) && return p

	return T(map(propertynames(p)) do k

		prop = getproperty(p,k)

		return k==:Weight ? prop*w : prop 

	end...)

end 

function Base.:*(p::T,w::Number)::T where T<:Union{Product,
																								 MixedProduct,
																								 Scalar}
	w*p 

end 



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function Base.:+(ps::Vararg{<:Product})::Scalar 

	+((MixedProduct(p) for p in ps)...)

end  


function Base.:-(p::P, q::Q)::Scalar where {T<:Union{Product,Scalar,MixedProduct},P<:T,Q<:T}

	p + (-1.0)*q 

end 



Base.sum(ps::AbstractVector{<:MixedProduct})::Scalar = Scalar(ps)

Base.:+(ps::Vararg{<:MixedProduct})::Scalar = Scalar(ps)


function Base.:+(S::Scalar, P::Product)::Scalar

	S + MixedProduct(P)

end  

function Base.:+(P::Product, S::Scalar)::Scalar

	MixedProduct(P) + S 

end  



function Base.:+(S::Scalar, P::MixedProduct)::Scalar

	S + Scalar(P)

end  

function Base.:+(P::MixedProduct, S::Scalar)::Scalar

	Scalar(P) + S

end 



function Base.:+(S1::Scalar, S2::Scalar)::Scalar


	new_Fields = union(fieldargs(S1),fieldargs(S2))


	W2 = [S2.Weight*weight(t2) for t2 in parts(S2)]

	nz2 = findall(!small_weight, W2)


	isempty(nz2) && return Scalar{length(new_Fields),length(S1)}(
																					Tuple(new_Fields), 
																					fieldargs_distrib(S1),
																					S1.Weight,
																					parts(S1)
																					)
									
	new_FieldNames = argnames.(new_Fields) 


	W1 = [S1.Weight*weight(t1) for t1 in parts(S1)]

	nz1 = findall(!small_weight, W1)


	if isempty(nz1) 
		
		new_FD = Tuple(fieldargs_distrib(new_FieldNames,argnames(S2,i)) for i=nz2)

		return Scalar{length(new_Fields),length(nz2)}(
																					Tuple(new_Fields),
																					new_FD,
																					S2.Weight,
																					Tuple(parts(S2,i) for i=nz2))
																					
	end 





	for i1 in nz1 

		for i2 in nz2 

			if !small_weight(W2[i2]) && proportional(parts(S1,i1), parts(S2,i2))

				W1[i1] += W2[i2]

				W2[i2] = 0   

				break 

			end

		end  

	end 


	nz1 = findall(!small_weight, W1)
	nz2 = findall(!small_weight, W2)


	new_FD = ((fieldargs_distrib(S1,i) for i=nz1)...,
						(fieldargs_distrib(new_FieldNames, argnames(S2,i)) for i=nz2)...)


	new_Terms = (
							(parts(S1,i)*(W1[i]/weight(S1,i)) for i=nz1)...,
							(parts(S2,i)*(W2[i]/weight(S2,i)) for i=nz2)...,
												)


	return Scalar{length(new_Fields),length(nz1)+length(nz2)}(
														Tuple(new_Fields), new_FD, 1.0, new_Terms)

end 
 


function Base.:*(S1::Scalar, S2::Scalar)::Scalar

	error("NOT yet checked")

	new_Fields = union(fieldargs(S1),fieldargs(S2))

	new_Terms = Tuple(p1*p2 for p1 in parts(S1) for p2 in parts(S2))

	#	new_FD  can be obtained directly from fieldargs_distrib(S1) and S2

	return Scalar_(new_Fields,
										fieldargs_distrib(new_Fields, new_Terms),
										S1.Weight*S2.Weight,
										new_Terms)

end 



















































































































































































































































#############################################################################
end 
