import myLibs: Algebra 

import QuadGK,Random


distr_names = [:Lorentzian,
							 :Gaussian,
							 :Rectangle,
							 :Heaviside]




area_under(f,w)::Float64 = round(QuadGK.quadgk(f, -50w,50w ,rtol=1e-14)[1],digits=2)


#@show w h 


verbose = false
N=100 


for F in distr_names
	
	
	verbose && println() 

	verbose && @info F 

	@testset "$F" begin 

	Random.seed!(Int(round(100*time())))

	for (w,h) in eachcol(rand(2,N).+0.3)

	
		verbose && println("\n ---------------------------------------------- ")


	
	
		verbose && println("\n ****** Only w is given")
	
		verbose && @show Algebra.get_bare_height_provided_width(F, w) Algebra.area_under_distrib(F, w)
	
	#	n1 = area_under(f1) 

#==========#
#		pf_w1 = 1/Algebra.area_under_distrib(F, w) 
		pf_w = Algebra.getNDistrib_prefactor(F,w)
#		@test pf_w1 ≈ pf_w
#==========# 

		f_w(x) = pf_w* Algebra.getDistrib(F, x, w)
	
		if F!=:Heaviside 
			
			verbose && @show area_under(f_w,w)
	
			@test isapprox(area_under(f_w,w),1,atol=F==:Lorentzian ? 0.03 : 1e-5)
	
		end 
	
		verbose &&	println("Max. height = ", f_w(F==:Heaviside ? w+1 : 0))
	
	
		verbose && println("\n ****** Both w and h are given")
		
#==========# 
#		pf_wh1 = h/Algebra.get_bare_height_provided_width(F, w)
		pf_wh = Algebra.getNDistrib_prefactor(F,w,h)
#		@test pf_wh1≈pf_wh
#==========# 

	
		f_wh(x) = pf_wh * Algebra.getDistrib(F, x, w) 
	
		verbose && (F==:Heaviside || @show area_under(f_wh,w))
	
		H_wh = f_wh(F==:Heaviside ? w+1 : 0)
	
		verbose && println("Max. height = ", H_wh)
	
		@test H_wh ≈ h
	
		
		
		
		verbose && println("\n ****** Only h=$h is given")
	
#==========# 
#		pf_h1, w_h1 = Algebra.get_width_provided_height(F, h) 
#		pf_h1 /= Algebra.area_under_distrib(F, w_h1)
		pf_h,w_h = Algebra.getNDistrib_prefactor_and_arg(F,h)
#		@test pf_h≈pf_h1 && w_h≈w_h1
#==========# 



		verbose && @show pf_h

		f_h(x) = pf_h * Algebra.getDistrib(F, x, w_h)


		if F!=:Heaviside 
			
			verbose && @show area_under(f_h,w_h)
	
			@test isapprox(area_under(f_h,w_h),1,atol=F==:Lorentzian ? 0.03 : 1e-5)
	
		end 

		H_h = f_h(F==:Heaviside ? w_h+1 : 0)
	
		verbose && println("Max. height = ", H_h)
	
		@test H_h ≈ h

	
		end 	


		println() 
	end 


end 


