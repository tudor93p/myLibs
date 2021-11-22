import myLibs:HoppingTerms ,H_Superconductor


HTs = [HoppingTerms.HoppingTerm(x,y) for x in [true,false] for y in [true, false]]


@testset "ordering" begin

for ht1 in HTs, ht2 in HTs 

	eq = ht1==ht2  

	l = ht1 < ht2 

	g = ht1 > ht2 

	ge = ht1 >= ht2 

	le = ht1 <= ht2 


	@test (eq|l)==le
	@test (eq|g)==ge


#	@test eq|l|g

	if ht1<ht2 && ht1.use_spin
		
		f = HoppingTerms.upgrade(ht1,ht2)


		g = H_Superconductor.ToBasis(true, zeros(2,2))

		x = rand(ComplexF64,2,2)

		@test isapprox(f(x),g(x))


	end 

end


end 
