module Algebra
#############################################################################



import ..LA

import ..Utils, ..ArrayOps, ..Groups

#import Dierckx,FFTW#,QuadGK 

import PolynomialRoots
 
const EPSILON = 1e-20

#julia> chi_improper(l,a) = (-1)^l * sin((l+1/2)a))/sin(a/2)
#julia> chi_improper(l,a) = (-1)^l * sin((l+1/2)a)/sin(a/2)

import ..Groups: PauliMatrix


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



function linearSystem_iter!(x::AbstractVector,
														M::AbstractMatrix,
														N::AbstractMatrix,
														b::AbstractVector,
														m::Int=1
														)::Nothing

	m==0 && return

	x .+= M\(N*x + b) 

	return linearSystem_iter!(x,M,N,b,m-1)

end 

														
function linearSystem_GaussSeidel!!(x::AbstractVector,
																 A::AbstractMatrix,
																 args...
																)::Nothing #AbstractVector

	x .= 0 

	N = -LA.triu(A,1)

	LA.tril!(A)

	linearSystem_iter!(x, A, N, args...)

end  
 

function linearSystem_GaussSeidel!(x::AbstractVector,
																 A::AbstractMatrix,
																 args...
																)::Nothing #AbstractVector

	x .= 0

	linearSystem_iter!(x, LA.LowerTriangular(A), -LA.triu(A,1), args...) 

end  

function linearSystem_GaussSeidel(A::AbstractMatrix{T1},
																	b::AbstractVector{T2},
																 args...
															 )::Vector{promote_type(T1,T2,Float64)
																				 } where {T1<:Number,T2<:Number}

	x = Vector{promote_type(T1,T2,Float64)}(undef, length(b))

	linearSystem_GaussSeidel!(x, A, b, args...)

	return x

end  

function linearSystem_Jacobi(A::AbstractMatrix{T1},
														 b::AbstractVector{T2},
															 args...
															 )::Vector{promote_type(T1,T2,Float64)
																				 } where {T1<:Number,T2<:Number}

	x = Vector{promote_type(T1,T2,Float64)}(undef, length(b))

	linearSystem_Jacobi!(x, A, b, args...)

	return x 
	
end  


function linearSystem_Jacobi!(x::AbstractVector,
															 A::AbstractMatrix,
															 args...
																)::Nothing #AbstractVector

	x .= 0

	M = LA.Diagonal(A)

	linearSystem_iter!(x, M, M-A, args...)

end  


function linearSystem_Jacobi!!(x::AbstractVector,
															 A::AbstractMatrix,
															 args...,
																)::Nothing #AbstractVector

	x .= 0


	M = copy(LA.Diagonal(A))

	A .= M-A 

	linearSystem_iter!(x, M, A, args...)

end  







#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function poly2roots_from3vals!(d::AbstractVector{Float64},
						c::AbstractVector{Float64},
					 x::Union{AbstractVector{<:Real},
										NTuple{3,Real}},
					 y::Union{AbstractVector{<:Real},
										NTuple{3,Real}},
					 refused_sol::Real,
					 )::Nothing


	circshift!(d, x, 1)

	circshift!(d, x - d, 1)
	
	@assert abs(prod(d))>1e-12

	d .*= y 

	c.= 0.0

	for (di,(j,k)) in zip(d,((2,3),(3,1),(1,2)))

		c[3] += di

		c[2] -= di * x[j]
		
		c[2] -= di * x[k]

		c[1] += di * x[j] * x[k] 

	end 

	d .= x

	for (i,root) in enumerate(PolynomialRoots.roots(c,view(d,1:2)))
														
		re = real(root) 

		d[i] = (re==root && x[1]<=re<=x[3]) ? re : refused_sol

	end 

	d[3] = refused_sol 

	return 

end 



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#



#
#function RotMat2!(R::AbstractMatrix{Float64}, theta::Real, Ax::Int=3; kwargs...)::Nothing 
#
#	@assert LA.checksquare(R)==2 
#
#	R[1,1] = cos(theta) 
#
#	R[2,2] = R[1,1]
#
#	R[1,2] = (-1)^Ax * sin(theta) 
#	
#	R[2,1] = -R[1,2]
#
#	return 
#
#end  

function RotMat2!(R::AbstractMatrix{Float64}, 
									theta::Real, 
									Ax::Int=3; 
								 kwargs...)::Nothing 

	Groups.O2Repr!(R, theta)

	if Ax==2 

		R[1,2] *= -1
		R[2,1] *= -1

	end 

	return 

end 



function RotMat(::Val{2}, theta::Real, args...; kwargs...)::Matrix{Float64}

	R = zeros(2,2) 

	RotMat2!(R, theta, args...)

	return R

end 


function RotMat(::Val{3}, theta::Real, Ax::Int; kwargs...)::Matrix{Float64}

	inds = filter!(!isequal(Ax), [1,2,3]) 

	R = zeros(3,3) 

	R[Ax,Ax] = 1.0

	RotMat2!(view(R,inds,inds), theta, Ax)

	return R

end



function RotMat(D::Int, theta::Real, Ax::Int=3; kwargs...)::Matrix{Float64}

	@assert 1<=Ax<=3 "Choose 1==x, 2==y, 3==z" 

	RotMat(Val(D), theta, Ax; kwargs...)

end 


function RotMat(theta::Real; kwargs...)::Matrix{Float64}

	RotMat(Val(2), theta; kwargs...)

end 



function RotVecs(vecs::AbstractMatrix{<:Real}, 
										rot::AbstractMatrix{Float64}; dim::Int
										)::Matrix{Float64}

	dim==2 && return rot*vecs 

	dim==1 && return vecs*transpose(rot)

	error()

end

#function RotVecs(v::AbstractVector{<:Real}, 
#										rot::AbstractMatrix{Float64}; dim::Int
#										)::Vector{Float64}
#
#	vec(RotVecs(Utils.VecAsMat(v, dim), rot; dim=dim))
#
#end 





function RotVecs(R::AbstractMatrix{Float64}; dim::Int, kwargs...)::Function 

	function rot_vecs(vecs::AbstractMatrix{<:Real})::Matrix{Float64}

		RotVecs(vecs, R; dim=dim)

	end

end 

function RotVecs(D::Int, args...; dim::Int, kwargs...)::Function 

	RotVecs(RotMat(D, args...; kwargs...); dim=dim)

end 


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function InvPartRatio(P::AbstractVector{<:Real})::Float64
	# ipr=1 means perfectly localized on one atom; 
	# ipr=length(P) means equally distributed on all sites

	S = sum(abs, P, init=0.0)
	
	return S<EPSILON ? length(P) : S^2/mapreduce(abs2, +, P)

end  

function InvPartRatio(Ps::AbstractArray{<:Real,N},
										 dim::Int)::Array{Float64,N-1} where N

	dropdims(mapslices(InvPartRatio, Ps, dims=dim), dims=dim)

end 






#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function Mean(f::Function, iter::Utils.List)

	mapreduce(f, +, iter)/length(iter) 


end 


function Mean(A::AbstractArray, dims::Union{Tuple,Int})::AbstractArray

	dims_ = unique(vcat(dims...))

	i = [d in dims_ ? 1 : Colon() for d in 1:ndims(A)]

	return sum(A, dims=dims)[i...]/prod(size(A)[dims_])

end


function Mean(A::Utils.List, dims::Nothing=nothing)

	sum(A)/length(A)

end




#===========================================================================#
#
# Linear independence of arrays
#
#---------------------------------------------------------------------------#


function Nr_LinearIndependentArrays(array_list)

	LA.rank(Array(hcat(reshape.(array_list,:)...)))

end

function LinearIndependentArrays(array_list)
	
	maxrank = Nr_LinearIndependentArrays(array_list)
	
	indep(inds) = Nr_LinearIndependentArrays(array_list[inds])==length(inds)


	function add_new(new_ind=2, inds=[1]) 

		(length(inds)==maxrank || new_ind==length(array_list)+1) && return inds 

		for i in inds

			!indep([i, new_ind]) && return add_new(new_ind+1, inds)
			
		end

		indep([inds;new_ind]) && return add_new(new_ind+1, [inds;new_ind])

		return add_new(new_ind+1, inds)

	end

	return matrix_list[add_new()]

end

#===========================================================================#
#
# Modified Gram-Schmidt
#
#---------------------------------------------------------------------------#

function GramSchmidt(V::T, iter_axis::Int=1;
										 normalize::Bool=true, tol::Float64=1e-10)::T where T

	inds(i) = insert!(repeat(Any[:], ndims(V)-1), iter_axis, i)

#	inds(i) = [j==iter_axis ? i : Colon() for j in 1:ndims(V)]


	iter = if T<:AbstractArray{<:Number}

							A->eachslice(A; dims=iter_axis)

					elseif T<:AbstractVector{<:AbstractArray{<:Number}}

							identity

					end 


	U = copy(V)*0.0

	for (j,Vj) in enumerate(iter(V))
		
		u = copy(Vj)

		for i in 1:j-1

			u -= (Ui->LA.dot(Ui,u)*Ui)(U[inds(i)...])

		end 

		norm = LA.norm(u)

		norm<tol && error("Invalid input matrix for Gram-Schmidt!")

		normalize && setindex!(U, u/norm, inds(j)...)

	end


  return U

end





#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function Commutator(A::AbstractMatrix, B::AbstractMatrix)::Matrix

	@assert size(A)==size(B) 

	A*B-B*A

end

function AntiCommutator(A::AbstractMatrix, B::AbstractMatrix)::Matrix

	@assert size(A)==size(B) 

	A*B+B*A

end



function Commutator(A::AbstractMatrix, B::AbstractMatrix,
										args::Vararg{Pair{Int,<:Function}}
										)::Matrix 

	@assert size(A)==size(B) 

	#a1,b1,b2,a2 = indexin(1:4, map(first,args))

	fs = Dict(args)

	return -(*(haskey(fs, 1) ? fs[1](A) : A,
						 haskey(fs, 2) ? fs[2](B) : B),
					 *(haskey(fs, 3) ? fs[3](B) : B,
						 haskey(fs, 4) ? fs[4](A) : A)
					 )


end




#===========================================================================#
#
# Orthonormal eigenvectors of a Hermitian matrix
#			grouped according to their energy (degenerate eigenvectors together)
#
#---------------------------------------------------------------------------#

function eigenvectors_degenerate_sectors(H; tol=1e-10)

	eigen = LA.eigen(LA.Hermitian(Array(H)))


	return [GramSchmidt(eigen.vectors[:,s],2) 
									for s in Utils.IdentifySectors(eigen.values; tol=tol)]

end

#===========================================================================#
#
#	Simultaneously diagonalize a set of commuting Hermitian matrices
#
#---------------------------------------------------------------------------#
	
project(Ms,U) = [U'*M*U for M in Ms]

function SimultDiagonaliz_CommutingHMatrices(Hs; tol=1e-10)

	all(length.(Hs).==1) && return Utils.UnitMatrix(1)

	projectors = eigenvectors_degenerate_sectors(Hs[1], tol=tol)


	length(Hs)==1 && return projectors



	return vcat(map(projectors) do p 
	
		np = SimultDiagonaliz_CommutingHMatrices(project(Hs[2:end], p); tol=tol) 

		return [p*new_p for new_p in np]

	end...)


end




#===========================================================================#
#
# Structure factors for a set of matrices
#
#---------------------------------------------------------------------------#

function structure_factors(Ls, test=true)::Array

	Utils.isList(Ls, AbstractMatrix) || error("List of matrices needed")

	F = zeros(Complex{Float64}, length(Ls), length(Ls), length(Ls))

	S = inv([LA.dot(L1,L2) for L1 in Ls, L2 in Ls])

	for (a,La) in enumerate(Ls)
		
		for (b,Lb) in enumerate(Ls[1:a-1])

			Cab = Commutator(La, Lb)

			F[a,b,:] = S*[-1im*LA.dot(L, Cab) for L in Ls]

			F[b,a,:] = -F[a,b,:]

		end

	end


	if test 

		e = map(eachcol(rand(axes(Ls,1),2,10))) do (a,b)
			
			LA.norm(Commutator(Ls[a],Ls[b]) .- im*mapreduce(prod,+,zip(F[a,b,:],Ls)))
		
		end

		maximum(e) > 1e-10 && println("####### Wrong structure factors")

	end

	return F

end






#===========================================================================#
#
# Computes the Killing form of a set of matrices 
# 													(or from the structure factors directly)
#
#---------------------------------------------------------------------------#

function Killing_form(arg::T) where T
	
	F = eachslice(if T<:AbstractArray{<:Number,3}

												arg
									
								elseif Utils.isList(T, AbstractMatrix) 

												structure_factors(arg)
								
								else 
												error("Wrong input for Killing form")
								end,

								dims=1)

	return [sum(Fa.*transpose(Fb)) for Fa in F, Fb in F]

end



#===========================================================================#
#
#	The center of a Lie algebra
#
#---------------------------------------------------------------------------#

function Center_LieAlgebra(Ls, K=Killing_form(Ls))

	Utils.isList(Ls, AbstractMatrix) || error("Wrong input") 

	return [mapreduce(prod,+,zip(l,Ls)) for l in eachcol(LA.nullspace(K))]

end


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function ConvoluteVectorPacket(weights, values, centers, 
															 delta, vectors::AbstractMatrix{<:Number};
															 dim::Int,
															 get_weights::Bool=false, normalize::Bool=true, 
															 keepdims::Bool=true, kwargs...) 

	W = getCombinedDistrib(weights, values, centers, delta;
												 normalize=normalize)

	out = dim |> function multipl(dim)
	
					v = sum(
									
							if dim==1 
						
						 W*view(vectors, axes(W,2),:)

							elseif dim==2 
								
								view(vectors, :, axes(W,2))*transpose(W)

							end,

							dims=dim)

							return keepdims ? v : dropdims(v, dims=dim)
					
					end 

	return get_weights ? (out,W) : out

end



#===========================================================================#
#
# Find root using bisection method
#
#---------------------------------------------------------------------------#

function FindRoot_Bisection(f,u10,u20=nothing;val10=nothing,val20=nothing,a=0,b=a+1e-8,move_right=x->x+1,move_left=x->x-1)

  function bisection(u1,u2)


 #   println("\nbisection")

 #   println("bisection ",u1," ",u2)

    um = (u1+u2)/2

    val = f(um)

    if val < a 
      return bisection(um,u2)

    elseif val > b
      return bisection(u1,um)

    else 
      return um

    end

  end



  uvalu(u) = (u,f(u))
 
# --- Find the interval, if not provided (properly) --- #

  function interval(u1,u2,val1,val2)

    (val1 < a) && (val2 > b) && return bisection(u1,u2)

#      println("bisection")
#      return bisection(u1,u2)
#    end

    if val1 > a

      if val1 > b
        (u2,val2) = (u1,val1)
      end

      (u1,val1) = uvalu(move_left(u1))
    end

    if val2 < b 

      if val2 < a
        (u1,val1) = (u2,val2)
      end

      (u2,val2) = uvalu(move_right(u2))
    end


#    println("interval ",u1," ",u2)

    return interval(u1,u2,val1,val2)

  end



  if isnothing(val10)

    val10 = f(u10)

  end




  if (isnothing(u20)  || isapprox(u10,u20,atol=1e-12))

    return interval(u10,u10,val10,val10)

  

  end

  return interval(u10,u20,val10,isnothing(val20) ? f(u20) : val20)

  


end





#===========================================================================#
#
# Partially invert matrix
#
#---------------------------------------------------------------------------#


function Cofactor(A::AbstractMatrix{T},row::Int,col::Int)::T where T

	(-1)^(col + row) * LA.det(view(A, axes(A,1).!=row, axes(A,2).!=col))

end



function Partial_Inverse(A,elements;compute="some")

  if compute == "some" 

    a,b = 1e-10,1e+10

    detA = LA.det(A)

    val = abs(detA)



    if a < val < b # && return [Cofactor(A,j,i) for (i,j) in elements]/detA

      return Dict([(i,j)=>Cofactor(A,j,i)/detA for (i,j) in elements])

    end


    Q1 = FindRoot_Bisection(u->abs(LA.det(A*u)),1,a=a,b=b,val10=val,
			move_right=x->x*2, move_left=x->x/2)

    Q2 = Q1^(1/(1-1/size(A,1)))

#    return A*Q2 |> B->[Cofactor(B,j,i) for (i,j) in elements]/LA.det(Q1*A)

    detA, B = LA.det(Q1*A), Q2*A

    return Dict([(i,j) => Cofactor(B,j,i)/detA  for (i,j) in elements])


#    return Partial_Inverse(A,elements,compute="full")


  elseif compute == "full" 

    return inv(Array(A)) |> B -> Dict([i=>B[i...] for i in elements])


  else
    
   error("'compute' must be 'full' or 'some'")

  end

end




#function Different_ijs(A,B=A)

#  A .!= reshape(B,1,A)
#
# 
#  for (i,a) in enumerate(A)
#
#    B[axes(B,1) .!= i]
# 
#
#  (j,
#
#  i,j = axes(a,1),axes(a,1)
#
#  i,j = repeat(i,inner=length(j)),repeat(j,outer=length(i))
#
#  w = (i.!=j)
# 
#
#

  
#  return [[y,x] for (x,y) in Base.product(B,A) if x!==y]




#end

#===========================================================================#
#
#  
#
#---------------------------------------------------------------------------#

function Cartesian_from_spherical(phi::Real)::Vector{Float64}

	Cartesian_from_spherical(1, phi)

end


function Cartesian_from_spherical(R::Real, phi::Real)::Vector{Float64}

	[R*cos(phi), R*sin(phi)]

end 

function Cartesian_from_spherical(R::Real, theta::Real, phi::Real
																 )::Vector{Float64}

	vcat(sin(theta)*Cartesian_from_spherical(R, phi), R*cos(theta))

end 


#===========================================================================#
#
# simulate localized wavefunctions 
#
#---------------------------------------------------------------------------#

#function Localized_WF(atoms,weight_orb,centers=atoms,localiz_shape="Gaussian",spread=1)
#
#  if ndims(centers) == 1 
#    centers = reshape(centers,1,:)
#  end
#
#  nr_orb = size(weight_orb(atoms[1,:]),1)
#
#  psi0 = zeros(Complex{Float64},size(atoms,1)*nr_orb,size(centers,1))
#
#  weights = get_Distribution(localiz_shape)
#
##  dist = OuterDist(atoms,centers)
#
#
#  W2 = weights.(eachcol(atoms),eachcol(centers),[spread])
#
#  W2 = ArrayOps.Normalize_Columns((*).(W2...))
##  W2[i,j] : weight on atom i for the wf with center j
#
#
#  for (iatom,(atom,w2)) in enumerate(zip(eachrow(atoms),eachrow(W2)))
#
#    w = ArrayOps.Normalize(weight_orb(atom)).*reshape(w2,1,:)
#
#    psi0[TBmodel.Hamilt_indices(1:nr_orb,iatom,nr_orb),:] = w
#
#  end
#
#  return psi0
#
#end

#===========================================================================#
#
# project a single wavefunction to certain eigenfunctions around E0 
#
#---------------------------------------------------------------------------#

function Project_WF_toE0(eigenen,eigenwf,newE0,newE0wid;projection_method="impose_weights",energy_shape="Gaussian")


  function change_coef()

    weights = get_Distribution(energy_shape)(eigenen,newE0,newE0wid)

    if projection_method == "impose_weights"

      return coef -> (coef./abs.(coef)) .* weights
    
    elseif projection_method == "partial_spectrum"
   
      weights2 = Complex{Float64}.(weights.>maximum(weights,dims=1)/100)

      return coef -> coef .* weights2
    
    end
    
  end


  newcoef = change_coef()

  return wf -> eigenwf*ArrayOps.Normalize_Columns(newcoef(eigenwf'*reshape(wf,:,1)))


end





function diff_vc(v::Tv, c::Tc) where {Tv,Tc}

	all(T->T<:AbstractArray,(Tv,Tc)) && return reshape(c,1,:) .- reshape(v,:,1)

	return c .- v

	error()

end 






#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

struct myDistrib

	F::Function 

	function myDistrib(f::Function; normalize=false) 
	
		hasmethod(f, (Real, Real)) || error("Wrong input function")
	
		g(x::Real, w::Real)::Real = f(x,w)
	
		g(v::Real, c::Real, w::Real)::Real = f(c-v, w)
	
		g(x::AbstractArray, w::Real)::Array = f.(x,w)
	
		g(v::AbstractArray, c::Real, w::Real)::Array = f.(c.-v, w)
	
		g(v::Real, c::AbstractArray, w::Real)::Array = f.(c.-v, w)
	
		g(v::AbstractVector, c::AbstractVector, w::Real)::Matrix = f.(Utils.VecAsMat(c,1) .- Utils.VecAsMat(v,2), w)
	


		function g(w::Real)::Function 
			
			h(x::Real)::Real = f(x,w)	
	
			h(v::Real, c::Real)::Real = f(c-v, w)
	
			h(x::AbstractArray)::Array = f.(x,w)
	
			h(v::AbstractArray, c::Real)::Array = f.(c.-v, w)	
	
			h(v::Real, c::AbstractArray)::Array = f.(c.-v, w)
	
			h(v::AbstractVector, c::AbstractVector)::Matrix = f.(Utils.VecAsMat(c,1) .- Utils.VecAsMat(v,2), w)
	
	
			return h 
	
		end 

	
		return new(g)

	
	end 
	
	
	function myDistrib(str::String; kwargs...)
	
		E = Meta.parse(str); 
		
		return myDistrib(@eval (x::Real, w::Real) -> $(E); kwargs...)
	
	end 


end 
# f = @eval ... ;return myDistrib(f)
#  Remembers only the last f and overwrites the previous ones
#
# myDistrib(@eval ...) might create strange "var is not defined" errors


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


normDistrib(W::Real)::Real = W 

function normDistrib(W::AbstractVector{T})::T where T<:Real 
	
	sum(W)

end 

function normDistrib(W::T)::T where T<:AbstractMatrix{<:Real} 

	sum(W, dims=2) # sums all centers => length(norm)==length(values)

end 

normDistrib(F::Function)::Function = normDistrib ∘ F 

normDistrib(D::myDistrib)::Function = normDistrib(D.F)


function normalizeDistrib(W::T, N::Float64
												 )::T where T<:Union{Real,AbstractVector{<:Real}}
	
	W/(N+EPSILON)

end 


function normalizeDistrib(W::AbstractArray{Tw,Nw}, norm::AbstractArray{Tn,Nn}
													)::AbstractArray{promote_type(Tw,Tn), max(Nw,Nn)} where {
																						Nw,Nn}  where {
																		Tw<:Number, Tn<:Number}

	W./(norm .+ EPSILON)

end 



function normalizeDistrib(W::Tw, N::Tn, n::Bool
												 ) where {Tw<:T,Tn<:T} where T<:Union{Real, AbstractVecOrMat}

	n ? normalizeDistrib(W, N) : W 

end 


function normalizeDistrib(W::T)::T where T<:Union{Real, AbstractVecOrMat}

	normalizeDistrib(W, normDistrib(W))

end 


function normalizeDistrib(W::T, n::Bool)::T where T<:Union{<:Real, <:AbstractVecOrMat}

	n ? normalizeDistrib(W) : W 

end 


normalizeDistrib(F::Function)::Function = normalizeDistrib ∘ F 
normalizeDistrib(F::Function, n::Bool)::Function = n ? normalizeDistrib(F) : F

normalizeDistrib(D::myDistrib)::Function = normalizeDistrib(D.F) 
normalizeDistrib(D::myDistrib, n::Bool)::Function = normalizeDistrib(D.F, n)





function normalizedDistribAndNorm(W::Union{Real, AbstractVecOrMat}, n::Bool=true)::Tuple

	N = normDistrib(W) 

	return (normalizeDistrib(W, N, n), N)

end 


function normalizedDistribAndNorm(F::Function, n::Bool=true)::Function 

	FN(args...) = normalizedDistribAndNorm(F(args...), n)

end 

function normalizedDistribAndNorm(D::myDistrib, n::Bool=true)::Function 

	normalizedDistribAndNorm(D.F, n)

end 





function (D::myDistrib)(args...; normalize::Bool=false) 

	normalizeDistrib(D.F(args...), normalize)

end 

#===========================================================================#
#
# Lorentzian, Gaussian functions etc
#
#---------------------------------------------------------------------------#


Lorentzian = myDistrib("w / (x^2 + w^2) / pi")

#Gaussian = myDistrib("exp(-x^2/w^2)/(w*sqrt(pi))")
Gaussian = myDistrib("exp(-0.5*x^2/w^2)/(w*sqrt(2*pi))")

Heaviside = myDistrib("Float64(x>=w)")

#Rectangle = myDistrib("Float64(-w/2<=x<=w/2)")
Rectangle = myDistrib("Float64(-w<=x<=w)")

FermiDirac = myDistrib("1/(exp(x/w)+1)")


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


#function get_width_provided_height(F::Symbol,h::Real)::NTuple{2,Float64}
#
#	F==:Lorentzian && return (1,1/(pi*h))
#
#	F==:Gaussian && return (1,1/(sqrt(pi)*h))
#
#	F==:Rectangle && return (1,1/h)
#
#	F==:Heaviside && return (h,0)
#
#end  
#
#
#function get_bare_height_provided_width(F::Symbol, w::Real)::Float64
#
#	F == :Lorentzian && return 1/(w*pi)
#
#	F == :Gaussian && return 1/(w*sqrt(pi))
#
#	return 1 
#
#end 
#
#
#function area_under_distrib(F::Symbol, w::Real)::Float64 
#
#	F==:Rectangle ? w : 1 
#	
#end 
#

function getNDistrib_prefactor(F::Symbol, w::Real,
															h::Nothing=nothing)::Float64 

	1/(F==:Rectangle ? 2w : 1)

end 

	
function getNDistrib_prefactor_and_arg(F::Symbol, h::Real
																	)::NTuple{2,Float64}

	pf, w = if F==:Lorentzian
																(1,1/(pi*h))
					elseif F==:Gaussian
																(1,1/(sqrt(2pi)*h))
					elseif F==:Rectangle
																(1,1/(2h))
					elseif F==:Heaviside
																(h,0)
					end 

	return pf*getNDistrib_prefactor(F, w), w
	
end 


function getNDistrib_prefactor(F::Symbol, w::Real, h::Real)::Float64

	h0 = 	if 			F==:Lorentzian	
																1/(w*pi)
				elseif 	F==:Gaussian 		
																1/(w*sqrt(2pi))
				else										
																1
				end 

	return h/h0
	
end 



#w=rand()+0.1;x=rand(100,100);
#@show LA.norm(Gaussian(x,w) - Lorentzian(x,w))
#
#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function getDistrib(args::Tuple; kwargs...)#::Function

	getDistrib(args...; kwargs...)

end 

function getDistrib(name::String, args...; kwargs...)#::Function
	
	getDistrib(Symbol(name), args...; kwargs...)

end 

function getDistrib(F::Function; normalize=false)::Function
	
	normalizeDistrib(F, normalize)

end

function getDistrib(name::Symbol, args...; kwargs...) 

	getDistrib(getproperty(@__MODULE__, name), args...; kwargs...)
	
end 


function getDistrib(F::Function, args...; normalize=false) 
	
	normalizeDistrib(F(args...), normalize)

end 

function getDistrib(D::myDistrib, args...; kwargs...)

	getDistrib(D.F, args...; kwargs...)

end 






function getDistribAndNorm(args...; kwargs...)

	normalizedDistribAndNorm(getDistrib(args...; normalize=false);
													 kwargs...)

end 






#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function getCombinedDistrib(args; normalize::Bool=false, kwargs...) 

	length(args)==1 && return getDistrib(args[1]...; normalize=normalize)

	D(a) = getDistrib(a..., normalize=false)

	D1 = D(args[1]) # may be a function or numeric


	if typeof(D1)<:Union{Real,AbstractArray}

		return normalizeDistrib(
							mapreduce(D, ArrayOps.multiply_elementwise, args[2:end], init=D1),
							normalize)


	elseif D1 isa Function 

		Ds = [D1; D.(args[2:end])]

		totalD(Args...) = mapreduce(d->d(Args...), ArrayOps.multiply_elementwise, Ds)

		return normalizeDistrib(totalD, normalize) 

	else 

		error("Wrong D1")

	end 

end 


function getCombinedDistribAndNorm(args; normalize::Bool=false)

	normalizedDistribAndNorm(getCombinedDistrib(args, normalize=false), 
													 normalize)

end 


#	If length('args')>1, the quantities must be given in the proper order. 
#	Either of:
#				names, deltas 
#				names, values, deltas 
#				names, values, centers, deltas 


mktuple(a) = isa(a,Tuple) ? a : tuple(a)

function getCombinedDistrib(args...; kwargs...)

	getCombinedDistrib(Utils.Zipmap(mktuple, args); kwargs...)

end


function getCombinedDistribAndNorm(args...; kwargs...)
	
	getCombinedDistribAndNorm(Utils.Zipmap(mktuple, args); kwargs...)

end 










#===========================================================================#
#
# solve Quadratic equation
#
#---------------------------------------------------------------------------#

function QuadraticEq(a::Number,b::Number,c::Number)::Vector

	PolynomialRoots.roots([c,b,a])

end



#===========================================================================#
#
# Cross product with vectors length 2
#
#---------------------------------------------------------------------------#

function cross(A::AbstractArray,B::AbstractArray)::AbstractArray


  good(r) = length(r) < 3 ? vcat(r,zeros(3-length(r))) : r[1:3]

  return LA.cross(good(A),good(B))


end


#===========================================================================#
#
# Dot product with vectors of different lengths
#
#---------------------------------------------------------------------------#

dot(A::Number, B::Number)::Number = dot([A],[B])
dot(A::Number, B::AbstractVector)::Number = dot([A],B)
dot(A::AbstractVector, B::Number)::Number = dot(A,[B])


function dot(A::AbstractArray{T,0}, B)::Number where T 
	
	dot(vcat(A...), B)

end 

function dot(A, B::AbstractArray{T,0})::Number where T 
	
	dot(A, vcat(B...))

end


function dot(A::AbstractVector, B::AbstractVector)::Number

  i = 1:min(length(A), length(B))

  return LA.dot(A[i], B[i])

end

#===========================================================================#
# 
# Efficient outer sum of two lists of vectors U = 
#	Returns a list [X,Y,...], such that U[i]+V[j] = X[i,j] + Y[i,j] + ...
#
#---------------------------------------------------------------------------#


function OuterBinary(U::Union{Number,AbstractArray}, op::Function;
										 kwargs...)::AbstractArray 

	OuterBinary(U, U, op; kwargs...)

end

function OuterBinary(U::Tu, V::Tv, op::Function; 
										 flat::Bool=false
										 )::AbstractArray where Tu<:T where Tv<:T where T<:Union{Number,AbstractVector}
	
	right_shape(a::Number) = hcat(a)
	right_shape(a::AbstractVector) = reshape(a,:,1)

	return dropdims(OuterBinary(right_shape(U), right_shape(V), op; 
															flat=flat, dim=1),
									dims=3-flat)

end 


function OuterBinary(U::Tu, V::Tv, op::Function; 
										 flat::Bool=false, dim::Int=1
										 ) where Tu<:T where Tv<:T where T<:AbstractMatrix

	dim2 = 3 - dim # columns if dims=rows and reversed


	sizes = map(enumerate(zip(size(U),size(V)))) do (d,(su,sv))

				if d==dim 

					return flat ? su*sv : (su,sv)

				else  
	
					su!=sv && error("Wrong sizes")

					return su 

				end 

			end 



	typ = Base.return_types(op, typeof.(first.((U,V)))) |> function aux(ts)

					!isempty(ts) ? ts[1] : typeof(op(U[1],V[1])) # Any

		end 



	out = similar(Array{typ}, Utils.flat(sizes)...)
	
	i0 = fill(Colon(), 2-flat)


	for i in 1:sizes[dim2]

		setindex!(out,
							view(op.(selectdim(U, dim2, i),
											 reshape(selectdim(V, dim2, i), 1, :)
											 ), i0...),
							[i0,i][dim]..., [i0,i][dim2]...)
		
	end 
	
	return out 

end


OuterSum(args...; kwargs...) = OuterBinary(args..., +; kwargs...)

OuterDiff(args...; kwargs...) = OuterBinary(args..., -; kwargs...)

function OuterDist2(args...; kwargs...) 
	
	D = get(kwargs, :dim, 1)==1 ? 3 : 1

	return selectdim(sum(abs2, OuterDiff(args...; kwargs...), dims=D),
									 D, 1) 
									
end  


function OuterDist(args...; kwargs...) 

	sqrt.(OuterDist2(args...; kwargs...))

end 


#===========================================================================#
# 
# Flattened outer sum, each line S[k] is a sum of U[i_k] and V[j_k]
#	S[k] = U[i] + V[j], with (i,j) = unravel_index(k,(len(U),len(V)))
#
#---------------------------------------------------------------------------#
  

FlatOuterSum(args...; kwargs...) = OuterBinary(args..., +; kwargs..., flat=true)

FlatOuterDiff(args...; kwargs...) = OuterBinary(args..., -; kwargs..., flat=true)

function FlatOuterDist(args...; kwargs...) 
	
	LA.norm.(eachslice(FlatOuterDiff(args...; kwargs...),
										 dims=get(kwargs,:dim,1)))

end 


#eachslice(D,dims=dim)
#
#
#def FlatOuter_IndexConvert(U,V):
#
#  get_k = lambda ij: np.ravel_multi_index(ij,(len(U),len(V)))
#
#  get_ij = lambda k: np.unravel_index(k,(len(U),len(V)))
#
#  return get_k,get_ij
#
#



#===========================================================================#
#
# Distance equals 
#
#---------------------------------------------------------------------------#

Dist(x::Real, y::Real)::Real = abs(x-y)

function Dist(x::Real, b::T)::T where {T<:AbstractVecOrMat}
	
	abs.(x-b)  

end 

function Dist(a::T, y::Real)::T where T<:AbstractVecOrMat
	
	abs.(a-y) 

end 

function Dist(a::AbstractVector{<:Real}, 
							b::AbstractVector{<:Real}; kwargs...)::Float64
	
	LA.norm(a-b)

end 


function Dist(a::AbstractVector{<:Real},
							B::AbstractMatrix{<:Real}; dim, kwargs...)::Vector{Float64}

	LA.norm.(eachslice(B.-reshape(a, insert!(Any[Colon()], dim, 1)...); 
										 dims=dim))

end 


function Dist(A::AbstractMatrix{<:Real},
							b::AbstractVector{<:Real}; kwargs...)::Vector{Float64}

	Dist(b, A; kwargs...)

end 



function Dist(A::AbstractMatrix{<:Real},
							B::AbstractMatrix{<:Real}; dim::Int, kwargs...)::Matrix{Float64}

	OuterDist(A, B; dim=dim)

end 


function Dist(t::Tuple{Vararg{<:Real}}, arg::T; kwargs...) where T<:Union{<:Real, AbstractVecOrMat{<:Real}}
	
	Dist(vcat(t...), arg; kwargs...)

end 

function Dist(arg::T, t::Tuple{Vararg{<:Real}}; kwargs...) where T<:Union{<:Real, AbstractVecOrMat{<:Real}}

	Dist(arg, vcat(t...); kwargs...)

end 

function Dist(t1::T1, t2::T2; kwargs...)::Float64 where {T1<:T,T2<:T} where T<:Tuple{Vararg{<:Real}}

	Dist(vcat(t1...), vcat(t2...); kwargs...)

end 



function EuclDistEquals(d0::Real; tol::Float64=1e-8, dim::Int,
											 kw...)::Function

	isd0(dist) = isapprox(dist, d0, atol=tol)

#	dim==1 means vectors on rows
# dim==2 means vectors on columns 




	return function EuclDistEquals_(A::a, B::b; kwargs...
																	) where {a<:T,b<:T} where T<:AbstractVecOrMat

		isd0.(Dist(A,B; dim=dim))

	end


end 


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


	
function LargestLengthscale(points::AbstractMatrix; dim::Int)::Float64

	Dist(Utils.zipmap(extrema, eachslice(points, dims=[2,1][dim]))...)

end 
	

#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#




function get_Bonds(atoms::AbstractMatrix, 
									 len_or_f::Union{<:Function,<:Real}; 
									 order_pairs::Bool=true,
									kwargs...)

	get_Bonds(atoms, atoms, len_or_f; order_pairs=order_pairs, kwargs...)
						
end 



function get_Bonds(atoms1::Union{Function,AbstractMatrix},
									 atoms2::Union{Function,AbstractMatrix}, 
									 bondlength::Real,
									 args...;
									 tol::Float64=1e-8, dim::Int, kwargs...)

	get_Bonds(atoms1, atoms2, 
						EuclDistEquals(bondlength; tol=tol, dim=dim),
						args...;
						dim=dim, kwargs...)

end 




#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function get_Bonds(get_atoms_batch::Function,
									 len_or_f::Union{Function,Real},
									 output_index::Function,
									 batches::AbstractVector{NTuple{2,Int}};
									 order_pairs::Bool=true,
									kwargs...)

	get_Bonds(get_atoms_batch, get_atoms_batch,
						len_or_f,
						output_index, output_index,
						batches;
						order_pairs=order_pairs,
						kwargs...)

end 



function get_Bonds(get_atoms_batch_1::Function,
									 get_atoms_batch_2::Function,
									 isBond::Function,
									 output_index_1::Function,
									 output_index_2::Function,
									 batches::AbstractVector{NTuple{2,Int}};
									 order_pairs::Bool=false,
									kwargs...)

	Utils.flatmap(batches) do (b,c)

		order = order_pairs && b==c 

		d = isBond(get_atoms_batch_1(b), get_atoms_batch_2(c); kwargs...)

		inds_batch_1, inds_batch_2 = output_index_1(b), output_index_2(c)

		return sort(map(findall(order ? LA.triu(d,1) : d)) do bond

			(i,j) = bond.I 

			I = inds_batch_1[i]
			J = inds_batch_2[j]

			return (!order || I<J) ? (I,J) : (J,I)
							
		end)

	end |> sort 

end 



function split_in_bathces(n::Int, N::Int)::Vector{OrdinalRange{Int,Int}}

	Utils.EqualDistributeBallsToBoxes_cumulRanges(n, cld(n,N))

end 





function get_Bonds(atoms1::AbstractMatrix, atoms2::AbstractMatrix, 
									 isBond::Function; 
									 dim::Int, order_pairs::Bool=false, N::Int=5000,
										)::Vector{NTuple{2,Int}}

	batches1 = split_in_bathces(size(atoms1,dim), N)
	batches2 = split_in_bathces(size(atoms2,dim), N)


	I = [(b,c) for b=axes(batches1,1)
			 for c=UnitRange(1, order_pairs ? b : length(batches2))]

	output_index_1(b::Int)::Vector{Int} = batches1[b]
	output_index_2(c::Int)::Vector{Int} = batches2[c]


	return get_Bonds(Utils.sel(atoms1, dim) ∘ output_index_1,
									 Utils.sel(atoms2, dim) ∘ output_index_2,
									 isBond,
									 output_index_1,
									 output_index_2,
									 I;
									 order_pairs=order_pairs,
									 dim=dim)

end 







#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

#
#sd_args1(a; dim::Int)::Tuple{AbstractVector{<:Real}, Int, Int} = (atoms1, dim, a)
#
#
#function bondRs_fromInds(bond_indices::AbstractVector{Tuple{Int,Int}}, 
#												 sd_args1::Function,
#												 sd_args2::Function=sd_args1
#												 )::Vector{Vector{Vector{<:Real}}}
#	f = Utils.sel(atoms1,dim) 
#
#	map(bond_indices) do (a,b)
#
#		collect(collect(f(i)) for (f,i) in zip((f1,f2),(a,b)))
#
#	end 
#end 

function bondRs_fromInds(inds_bonds::AbstractVector{NTuple{2,Int}},
												 get_atoms_batch::Function,
												 output_index::Function,
												 batches::AbstractVector{NTuple{2,Int}}; 
												 same_batches::Bool=true,
												 kwargs...
												 )::Vector{Vector{Vector{Real}}}

	bondRs_fromInds(inds_bonds,
									get_atoms_batch, get_atoms_batch,
									output_index, output_index,
									batches;
									same_batches=same_batches,
									kwargs...)

end 


bR_get_ind(a::Int, ::Colon)::Tuple{Bool,Int} = (true,a)  
bR_get_ind(::Int, ::Nothing)::Tuple{Bool,Int} = (false,0) 
bR_get_ind(::Int, A::Int)::Tuple{Bool,Int} = (true,A)
bR_get_ind(a::Int,v::AbstractVector{Int})::Tuple{Bool,Int} = bR_get_ind(a, only(indexin(a, v)))

function bondRs_fromInds(inds_bonds::AbstractVector{NTuple{2,Int}},
												 get_atoms_batch_1::Function,
												 get_atoms_batch_2::Function,
												 output_index_1::Function,
												 output_index_2::Function, 
												 batches::AbstractVector{NTuple{2,Int}}; 
												 same_batches::Bool=false,
												 kwargs...
												 )::Vector{Vector{Vector{Real}}}

	bondRs = [Vector{Vector{Real}}(undef,2) for bond in inds_bonds] 


	current = BitVector(undef, length(inds_bonds))

	As = zeros(Int, length(inds_bonds))

	Bs = zeros(Int, length(inds_bonds))



	for sector in Utils.IdentifySectors(first, batches)

		sub_batch = batches[sector] 

		i1 = sub_batch[1][1]

		atoms1 = get_atoms_batch_1(i1) 

		inds_batch_1 = output_index_1(i1) 


		for (i1,i2) in sub_batch 

			for (bond,((a,),Rs)) in enumerate(zip(inds_bonds,bondRs))

				if isassigned(Rs)

					current[bond] = false 

				else 

					current[bond], As[bond] = bR_get_ind(a, inds_batch_1) 

				end 

			end 


			potential_bonds = findall(current) 

			isempty(potential_bonds) && continue  



			same = same_batches && i1==i2 

			inds_batch_2 = same ? inds_batch_1 : output_index_2(i2)  


			for bond in potential_bonds 

				current[bond], Bs[bond] = bR_get_ind(inds_bonds[bond][2], 
																						 inds_batch_2) 
			end 


			potential_bonds = potential_bonds[current[potential_bonds]]

			isempty(potential_bonds) && continue  


			atoms2 = same ? atoms1 : get_atoms_batch_2(i2)


			for bond in potential_bonds 

				bondRs[bond] .= bondRs_fromInds((As[bond],Bs[bond]), atoms1, atoms2;
																				kwargs...)

			end 
	
		end # sub-batch

			#for ((a,b), Rab) in zip(inds_bonds, bondRs)
#	
#				isassigned(Rab) && continue 
#
#				A,success = bR_get_ind(a, inds_batch_1)
#
#				success || continue 
# 
#				B,succes = bR_get_ind(b, inds_batch_2)
#				
#				success || continue  
#
#				Rab .= bondRs_fromInds((A,B), atoms1, atoms2; kwargs...) 
#
#			end 
#
#		end 

	end # sector 


	@assert all(isassigned, bondRs) 
	
	return bondRs 


end 

	#Utils.flatmap(batches) do (b,c)

	#	order = order_pairs && b==c 

	#	d = isBond(get_atoms_batch_1(b), get_atoms_batch_2(c); kwargs...)

	#	inds_batch_1, inds_batch_2 = output_index_1(b), output_index_2(c)

	#	return sort(map(findall(order ? LA.triu(d,1) : d)) do bond

	#		(i,j) = bond.I 

	#		I = inds_batch_1[i]
	#		J = inds_batch_2[j]

	#		return (!order || I<J) ? (I,J) : (J,I)
	#						
	#	end)

	#end |> sort 





function bondRs_fromInds(bond_indices::AbstractVector{Tuple{Int,Int}}, 
												 atoms::AbstractMatrix{<:Real};
												 same_batches::Bool=true,
												 kwargs...)::Vector{Vector{Vector{<:Real}}}

	bondRs_fromInds(bond_indices, atoms, atoms; 
									same_batches=same_batches,
									kwargs...)

end 

function bondRs_fromInds(bond_indices::AbstractVector{Tuple{Int,Int}}, 
												 atoms1::AbstractMatrix{<:Real},
												 atoms2::AbstractMatrix{<:Real};
												 kwargs...)::Vector{Vector{Vector{<:Real}}}
												 
	[bondRs_fromInds(I, atoms1, atoms2; kwargs...) for I in bond_indices]

end 



function bondRs_fromInds(I::NTuple{2,Int},
												 atoms::Vararg{AbstractMatrix{<:Real},2};
												 dim::Int,
												 kwargs...)::Vector{Vector{<:Real}}

	[collect(selectdim(adi...)) for adi in zip(atoms, (dim,dim), I)]

end 




#function output_get_Bonds(bond_indices, atoms1, atoms2=atoms1; 
#													dim=1, inds=true, pos=false, as_matrices=false)

#	as_matrices && return get_Bonds_toMatrix(atoms1, bond_indices, atoms2;
#																					 inds=inds, pos=pos, dim=dim)
#
#	!pos && return bond_indices
#
#	bond_Rs = bondRs_fromIndices(bond_indices, atoms1, atoms2; dim=dim)
#
#	return !inds ? bond_Rs : (bond_indices, bond_Rs)
#
#	#	Rs[pair_index][member_index][coordinate_index]

#end


function get_Bonds_toMatrix(X; inds=false, pos=false, dim=1)::Matrix

	xor(inds,pos) || error("Specify either 'inds' or 'pos'")

	M = ArrayOps.init_zeros(X, dim=>length(X), [2,1][dim]=>sum(length, X[1]))

	for (i,x) in enumerate(X)

		setindex!(selectdim(M, dim, i), vcat(x...), :)

	end 

	return M

end 


function get_Bonds_asMatrix(args...; inds=false, pos=false, kwargs...)::Matrix

	get_Bonds_toMatrix(get_Bonds(args...; kwargs...), inds=inds, pos=pos)

end 


function get_Bonds_toMatrix(atoms1::AbstractMatrix, 
														bond_indices::Vector{Tuple{Int,Int}},
														atoms2::AbstractMatrix=atoms1; 
														dim=1, kwargs...)::Matrix

	get_Bonds_toMatrix(atoms1,
										 get_Bonds_toMatrix(bond_indices; inds=true,dim=dim),
										 atoms2;
										 dim=dim,kwargs...)

end 

function get_Bonds_toMatrix(atoms1::AbstractMatrix, 
														M_bond_indices::AbstractMatrix{Int},
														atoms2::AbstractMatrix=atoms1; 
														inds=false, pos=false, dim=1)::Matrix

	xor(inds,pos) || error("Specify either 'inds' or 'pos'")

	!pos && inds && return M_bond_indices


	dim2 = [2,1][dim]

	return cat(map(enumerate((atoms1,atoms2,))) do (i,atoms)

		selectdim(atoms, dim, selectdim(M_bond_indices, dim2, i)) 

	end..., dims=dim2)


end 




function get_Bonds_fromMatrix(M; inds=false, pos=false, dim=1)

	xor(inds,pos) || error("Specify either 'inds' or 'pos'")
	
	inds && return Tuple.(eachslice(convert(Array{Int},M),dims=dim))

	pos && return map(eachslice(M,dims=dim)) do Rij 

						n = div(length(Rij),2)

						return [Rij[1:n],Rij[n+1:end]]
				
				end 

end 	


#############################################################################

end
