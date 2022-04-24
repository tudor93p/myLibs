import myLibs: Utils 
#import Device: GL 
import PyPlot,LinearAlgebra
using Constants: VECTOR_STORE_DIM 
import QuadGK ,Random,Combinatorics

import Device: utils, Taylor 


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function rand_weight()::Union{Float64,ComplexF64} 
	
	rand([Float64,ComplexF64])(rand(0:10)/10)

end 

function rand_field(p) 

	rand(rand([Float64,ComplexF64]),
			 fill(3,Taylor.field_rank(p))...)
	
end


function test_derivative(original, deriv, field::AbstractArray, p, k::Int, I::Taylor.Index, allfields)

	z(x) = [k1==k ? x : x0 for (k1,x0) in enumerate(allfields)]

	F(x) = original(z(x)...)

	D(x) = deriv(z(x)...)

	expstep(dx) = onehot(p, dx, I)

	return utils.test_derivative(F, D, field, expstep)

end 

many_inds = [ [ 1, [1], (1,), [(1,)], ([1],), ],
							([1,2],(1,2),[(1,),(2,)],([1],[2])),
						 ([1,3],),
								(	[(1,2)],([1,2],) ),
								[ [(1,2,3),(2,3,1)], ([1,2,3],[2,3,1]) ],
								[ [(2,1,3),(2,1,3)], ],
								]


many_nr_rank = [(1,1),(2,1), (2,1), (1,2),(2,3), (2,3)
								]

#onehot(n::Int, i::Int)::Vector{Float64} = onehot(n, 1.0, i)  
#
#function onehot(n::Int, val::T, i::Int)::Vector{T} where T<:Number  
#
#	setindex!(zeros(T,n), val, i)
#
#end 

function onehot(p, dx::Float64, q_::Taylor.Index)

	setindex!(zeros(ComplexF64, fill(3,Taylor.field_rank(p))...), dx, q_.I...) 

end 
#flat_prod(ns) = collect(Base.Iterators.flatten(Base.product(ns...)))
flat_prod(ns) = collect(Base.product(ns...))[:]


rand_items(it, m::Int) = Utils.Random_Items(it, min(m,length(it)))  

flat_prod(ns, m::Int) = rand_items(flat_prod(ns), m)

sorted_prod(args...) = filter(issorted, flat_prod(args...))

function combs(n::Int)#,m::Int)

#	Random.seed!(rand(100:500)+Int(10round(time()))+ n)

#	filter(issorted, flat_prod(fill(axes(many_inds,1), n)))
	sorted_prod(fill(axes(many_inds,1), n))

end 

combs(n::Int, m::Int) = rand_items(combs(n), m)
														 

@time @testset "struct product basic" begin   

	@test !Taylor.proportional(Taylor.Product(1,1),Taylor.Product(1,2))


	
	for (S,inds_) in zip(many_nr_rank,many_inds)


#			Random.seed!(rand(100:500)+Int(10round(time()))+ rand(S)+sum(sum,rand(inds_)))

			M1 = Taylor.parse_inds(inds_[1]...)
	
			@test (Taylor.nr_fields(M1),Taylor.field_rank(M1))==S  
	
			P1 = Taylor.Product(rand_weight(),M1)
			
			@test zero(P1)(rand_field(P1))≈0



			!Taylor.small_weight(P1)
		
			for inds in inds_
			
				M2 = Taylor.parse_inds(inds...)
				
				@test (Taylor.nr_fields(M2),Taylor.field_rank(M2))==S  




				@test M2==M1
		
				@test Taylor.proportional(P1, Taylor.Product(rand_weight(),inds...))
			
				@assert Taylor.proportional(P1, Taylor.Product(rand_weight(),inds...))
				
		
			end  

			@test length(unique(Taylor.each_fieldfactor(P1)))==length(Taylor.count_unique(P1))
				
	end 

end 

println("\n"); #error()




@time @testset "struct product derivative" begin  

	for (S,inds_) in zip(many_nr_rank,many_inds) 

#		break
	
#			Random.seed!(rand(100:500)+Int(10round(time()))+ rand(S)+sum(sum,rand(inds_)))

		M1 = Taylor.parse_inds(inds_[1]...)

		P1 = Taylor.Product(rand_weight(),M1)
		
		ps1 = [Taylor.Product(1.0*i, I...) for (i,I) in enumerate(inds_)]

#		@show length(ps1) 
		ps2 = Taylor.cumulate(ps1)

		ps3 = Taylor.cumulate_(ps1) 

#		@show length(ps2) 
#		 Taylor.cumulate(ps1)
		
		@test utils.has_disjoint_pairs(==, ps2, ps3)


#		@show Taylor.cumulate_categories(ps1)
#		@time Taylor.cumulate_categories(ps1)

#		@show length(Taylor.cumulate_(ps1))
#		@time length(Taylor.cumulate_(ps1));

		@test length(ps2)==1 
		
		@test Taylor.proportional(ps2[1], P1)
	
		@test only(ps2).Weight == div(length(inds_)*(length(inds_)+1),2)



#		println("\n*** Term: ",join(inds_[1],","),"  ",P1.Weight,"  nr_fields=",Taylor.nr_fields(P1))

		J1,Q1 = Taylor.derivatives(P1)

		for (j1,q1) in zip(J1,Q1)

		
			@test Taylor.derivative(P1, j1)==q1

		end 




		cu = Taylor.count_unique(P1)


		for q in Base.product(fill(1:3, Taylor.field_rank(P1))...)

			q_ = Taylor.Index(q)
#			q_ = Tuple(q...)


			d = Taylor.derivative(P1, q_)
	
			@test utils.test_derivative(P1, d, rand_field(P1), dx -> onehot(P1, dx, q_)) 

			if !any(==(q_),M1)
				
#				@show P1 
#				@show q_ 
#				@show d 

				@test Taylor.small_weight(d)

			else 

#				@show q_  d

				@test Taylor.nr_fields(d)==Taylor.nr_fields(P1)-1


				j1 = findall([c==q_ for c in J1])

				@test length(j1)==1

				@test cu[only(j1)][2]*P1.Weight≈d.Weight



					

				if Taylor.nr_fields(d)==0 

					field = rand(ComplexF64, fill(3,Taylor.field_rank(P1))...)

					@test d(field)≈d.Weight 

				end

#				println()
			end 

		end

#	end  
	
#	for trial in 1:3

		#break 
#		Random.seed!(trial+Int(10round(time())))
	
		ps = vcat([[Taylor.Product(10.0^i, inds...) for inds in inds_[Utils.Random_Items(1:length(inds_))]] for (i,inds_) in enumerate(many_inds)]...)

		ps_ = Taylor.cumulate(Vector{Taylor.Product}(ps[Random.shuffle(1:length(ps))]))

		@test length(ps_)==length(many_inds)

		@test Set([Int(floor(log10(p.Weight))) for p in ps_]) == Set(1:length(many_inds))

	end 
end 

println("\n"); #error()


@time @testset " mixed product struct (old) " begin 

	for ((n,r),inds) in zip(many_nr_rank,first.(many_inds))


		P1 = Taylor.Product("eta",rand_weight(),inds...)

		for ((ncc,rcc),indscc) in zip(many_nr_rank,first.(many_inds))

			r==rcc || continue   

			P2 = Taylor.Product("eta*",rand_weight(),indscc...)


			P = Taylor.MixedProduct(("eta","eta*"), P1, P2)


			@test zero(P)(rand_field(P1),rand_field(P2))≈0


			#(P1.Weight,inds...),(P2.Weight,indscc...))

			for X in (2P,3P, 4P*P1, 5P1*P, 6P2*P,7P*P2)

				@test Taylor.argnames(X) == ("eta","eta*")
			
				@test only.(Taylor.argnames.((X,),(1,2))) == ("eta","eta*")

			end 

			@test only(unique(Taylor.field_rank(P)))==Taylor.field_rank(P1)==Taylor.field_rank(P2)

			@test sum(Taylor.nr_fields(P))==Taylor.nr_fields(P1)+Taylor.nr_fields(P2)

			field = rand(ComplexF64, fill(3,Taylor.field_rank(P1))...)
			fieldcc = rand(ComplexF64, fill(3,Taylor.field_rank(P2))...)

			@test P(field,fieldcc) isa Number

			P(field,fieldcc) == P.Weight*P1(field)*P2(fieldcc)

			for i=1:Taylor.nr_fields(P1)

				Taylor.derivative(P,1,i)

			end  

			for j=1:Taylor.nr_fields(P2)

				Taylor.derivative(P,2,j)

			end 

			for (field,field2) in [(1,2),(2,1)]#in [:field, :fieldcc]

				for q in Base.product(fill(1:3, Taylor.field_rank(P,field))...)

					q_ = Taylor.Index(q)
	
					d = Taylor.derivative(P, field, q_)  
			
					@test Taylor.argnames(d,1)==("eta",)
					@test Taylor.argnames(d,2)==("eta*",)


					if q_ in P.Factors[field] 

						@test Taylor.nr_fields(d,field)+1==Taylor.nr_fields(P,field)==[n,ncc][field]

						@test Taylor.nr_fields(d,field2)==Taylor.nr_fields(P,field2)==[n,ncc][field2]
					else 

						@test Taylor.small_weight(d)

					end 
	
				end 

			end 	

		end 

	end 

end 






println("\n"); #error()



@time @testset "struct mixed product basic + deriv" begin  

	 @test !Taylor.proportional(Taylor.MixedProduct(Taylor.Product(1,1)),
												 Taylor.MixedProduct(Taylor.Product(1,2)))


	 for n=[1,3]

		 #break 

		for iii in combs(n, 2)

			ns = Utils.Random_Items('a':'z',length(iii))

			prods = [Taylor.Product(fn, rand_weight(), rand(many_inds[i])...) for (fn,i) in zip(ns,iii)]

			w = rand_weight()



			@test 			w isa Union{Float64,ComplexF64}

			@test prods isa Union{AbstractVector{<:Taylor.Product},
																	Tuple{Vararg{Taylor.Product}}}

			@test Tuple(prods) isa Union{AbstractVector{<:Taylor.Product},
																	Tuple{Vararg{Taylor.Product}}}
		
			 P = Taylor.MixedProduct(w, prods) 


#			!Taylor.small_weight(P)

			fields = rand_field.(prods)

			a1 = prod([p(f) for (p,f) in zip(prods,fields)])
			
		 	A2 = P(fields...)

				A1 = w*a1 

			@test Taylor.same_weight(A1,A2)
			
			 a3 = prod(prods)(fields...)

			@test a1≈a3 
		
			A4 = (w*prod(prods))(fields...)

			@test Taylor.same_weight(A1,A4)

			A5 = prod(vcat(prods[1]*w,prods[2:end]))(fields...)

			@test Taylor.same_weight(A1,A5)


			i0 = rand(1:n) 

 			P2 = Taylor.MixedProduct_sameFields(P, prods[i0], i0) 

		 	@test Taylor.same_weight(P2(fields...),A1* prods[i0](fields[i0]))



#			@show P 
			for deriv_field in 1:n 
				
				J1,Q1 = Taylor.derivatives(P,deriv_field) 

				p7 = P.Factors[deriv_field] 

				for (factor,q_) in Taylor.enumerate_fieldfactors(p7)

#					D2 = only(Q1[[q_==j1 for j1 in Taylor.each_fieldfactor(J1)]])

		 			D = Taylor.derivative(P, deriv_field, q_)#factor) degeneracy!

					f3 = Taylor.fieldargs(P)[deriv_field]

					@test D==Taylor.derivative(P, f3, q_)
		
					fake_field = Taylor.FieldPlaceholder{25}("adasdsad")
					fake_inds = Taylor.Index(ntuple(i->rand(1:3),25))

					@test zero(P)==Taylor.derivative(P, fake_field, fake_inds)
					
					fake_inds = Taylor.Index(ntuple(i->rand(100:200),Taylor.field_rank(f3)))
					
					@test zero(P)==Taylor.derivative(P, f3, fake_inds)




					@test test_derivative(P,D,rand_field(p7),p7,deriv_field,q_,
																rand_field.(P.Factors))

				end 


				
				#continue 

				for (j1,q1) in zip(J1,Q1)

			
					X = [Taylor.derivative(P, deriv_field, i) for i in findall(j1, P.Factors[deriv_field])] 
					
					D1 = Taylor.cumulate(X)

					if length(D1) == 1 

						@test only(D1)==q1 

					else 


						@test all(Taylor.small_weight,X)

					end 


				end 

			end 


	#for (S,inds_) in zip(many_nr_rank,many_inds)
	
#			Random.seed!(rand(100:500)+Int(10round(time()))+ rand(S)+sum(sum,rand(inds_)))


		end 

	end 
end 





println("\n");	




function pr(n::AbstractString,s2)

	println("\n*** $n ***")

	for k in propertynames(s2)

		println("$k: ",getproperty(s2,k))

	end 

	applicable(Taylor.weight, s2) && 	@show Taylor.weight(s2)


	println() 

end 




@time @testset "scalar struct basics" begin 

#	Random.seed!(1)

	for nr_terms in [1,3], term_length in flat_prod(fill(1:3,nr_terms),3)

		#break 
#		@show term_length
		for iiiii in flat_prod(combs.(term_length,2), 2)

			possib_names = [Utils.Random_Items(k,rand(1:length(k))) for k in [
												["A","B"],
											["CC","DD"],
											["EEE","FFF"],
											]]


			prods = map(iiiii) do iii

				map(iii) do i 

					w = rand_weight() 
					
					I = rand(many_inds[i])
					
					fn = rand(possib_names[Taylor.field_rank(Taylor.Product(I...))])

					return Taylor.Product(fn, w, I...)

				end 

			end 

			@test nr_terms==length(prods)==length(iiiii)

			for (L1,L2) in zip(length.(prods),term_length)

				@test L1==L2 

			end 


			terms = [Taylor.MixedProduct(rand_weight(), vcat(p...)) for p in prods]


			for ps in prods 

				t1 = Taylor.MixedProduct(rand_weight(), vcat(ps...)) 
				
				t2 = Taylor.MixedProduct(rand_weight(), reverse(ps)...)
			
				if !all(Taylor.proportional(z1,z2) for (z1,z2) in zip(t1.Factors, t2.Factors)) 

					@test Taylor.proportional(t1,t2)

				end 
	
			end 
		
			
			S = Taylor.Scalar(rand_weight(), terms) 


			fields = Dict()

			for (rk,pn) in enumerate(possib_names)

				for n in pn 

					fields[n] = rand(ComplexF64, fill(3,rk)...)

				end 

			end 

			@test zero(S)([fields[k] for k in Taylor.argnames(S)]...)≈0


			z = Taylor.Scalar_(Taylor.fieldargs(S),
								 Taylor.fieldargs_distrib(S),
								 0.0,
								 Taylor.parts(S)
								 )

			@test z([fields[k] for k in Taylor.argnames(S)]...)≈0 

			@test zero(S)==z



			#continue
			
			@test S==Taylor.Scalar(S.Weight, terms...)

			@test allunique(Taylor.argnames(S))

#			@show allunique(mapreduce(Taylor.argnames, vcat, terms))

			for (i,t) in zip(S.FieldDistrib,S.Terms)

				@test Taylor.argnames(S)[i]==Taylor.argnames(t)

			end 


			reps = sum(1:length(S.Terms)-1; init=0) do i 

				sum(i+1:length(S.Terms); init=0) do j 

					length(intersect(Taylor.argnames(S, i),Taylor.argnames(S, j)))

				end 

			end   



			
			val1 = S((fields[k] for k in Taylor.argnames(S))...)

			val2 = sum(zip(prods,terms)) do (ps,term)

				P = Taylor.MixedProduct(term.Weight, (ps[i] for i in Random.shuffle(reverse(1:length(ps))))...)

				@test P==term 


				ft = (fields[k] for k in Taylor.argnames(term))


				@test term(ft...)≈P((fields[k] for k in Taylor.argnames(P))...)

				return (term*S.Weight)(ft...)

			end 


			@test Taylor.same_weight(val1, val2)

			s1 = +(terms...) 
			
			@test s1 isa Taylor.Scalar
		
			s2 = sum(terms) 
			
			@test s2 isa Taylor.Scalar

			@test s2==s1 

			s3 = sum(S.Terms)

			@test s3 isa Taylor.Scalar 

			Taylor.small_weight(S.Weight) || @test Taylor.equal_prods(s2,s3)
			
			
			if !Taylor.small_weight(S.Weight) && !Taylor.equal_prods(s2,s3)



					println("\n --- Terms --- ")
	
				pr.("t".*string.(1:length(terms)),terms)

				pr("s1",s1)
				pr("s2",s2)
				pr("s3",s3)
			
				pr("S",S)
			pr.("S Term ".*string.(1:length(S)),S.Terms)

				error()

			end 




		
			S3  = s3*S.Weight 


			@test Taylor.equal_prods(S,S3)
			

			S4 = S+S 

			S5 = 2S 

			@test S4==S5  


			@test S + 2S + 0.5S == S*3.5

		end 

	end 

end 

println("\n")

@time @testset "scalar derivative and tensors" begin 

#	Random.seed!(1)

	for nr_terms in 1:3, term_length in flat_prod(fill(1:3,nr_terms),3)

		#break 

		for iiiii in flat_prod(combs.(term_length,2), 2)

			possib_names = [Utils.Random_Items(k,
																				 rand(1:length(k)),
																				 ) for k in [
												["A","B"],
											["CC","DD"],
											["EEE","FFF"],
											]]


			prods = map(iiiii) do iii

				map(iii) do i 

					w = rand_weight() 
					
					I = rand(many_inds[i])
					
					fn = rand(possib_names[Taylor.field_rank(Taylor.Product(I...))])

					return Taylor.Product(fn, w, I...)

				end 

			end 


			

			for p in prods 

				args = Taylor.argnames(Taylor.MixedProduct(rand_weight(), p))

#				allunique(Taylor.argnames.(p))  

				@test allunique(args)
			
				@test args == Taylor.argnames(Taylor.MixedProduct(p))
				
				@test args == Taylor.argnames(Taylor.MixedProduct(p...))
				
				@test args == Taylor.argnames(Taylor.MixedProduct(rand_weight(), p...))
				
			#	@show Taylor.argnames(Taylor.MixedProduct_(rand_weight(), p))


			end  

			#mprods = [Taylor.MixedProduct(rand_weight(), p) for p in prods]


			S = Taylor.Scalar(rand_weight(), [Taylor.MixedProduct(rand_weight(), p) for p in prods])


			fields = Dict()

			for (rk,pn) in enumerate(possib_names)

				for n in pn 

					fields[n] = rand(ComplexF64, fill(3,rk)...)

				end 

			end  



			fs = [fields[k] for k in Taylor.argnames(S)]

			@test S(fs...)*0 ≈ zero(S)(fs...)

			for (k,f) in enumerate(Taylor.fieldargs(S))

				nzI = Taylor.Index.(unique(vcat([Taylor.each_fieldfactor(Taylor.parts(P,i4)) for P in Taylor.parts(S) for i4 in findall(==(f), Taylor.fieldargs(P))]...)))

				pzI = Taylor.Index.(eachcol(rand(1:3, Taylor.field_rank(f), 10)))

				for I in unique(vcat(pzI..., nzI...))

					@test I isa Taylor.Index
					
					D = Taylor.derivative(S, f, I)
	
					if !in(I,nzI) 
						
						@test D==zero(S)

						@test length(D)==1 
					end 

					@test xor(!(S(fs...)≈0) && in(I,nzI), D(fs...)≈0)

					@test test_derivative(S, D, fs[k], f, k, I, fs)

				end 


				inds,scalars = Taylor.derivatives(S,f)

				if isempty(inds) 
					inds = [Taylor.Index(fill(1,Taylor.field_rank(f)))]
					scalars = [zero(S)]
				end 

				T = Taylor.Tensor_fromDeriv(S,f)

				@test ndims(T)==1 

				@test size(T,1)==Taylor.field_rank(f)

				for (i,s) in zip(inds,scalars)

					@test T[i](fs...)≈ s(fs...)

				end 

				for (((i1,),t1),i2) in zip(T,inds)

					@test i1==i2

				end 


			end 




		end 

	end 

end 






#

	
	
	
	
	
	
