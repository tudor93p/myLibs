using Revise,Test,Distributed 
using SharedArrays
import myLibs:Utils,MeshInPlace

NR_ENERGIES= 400


@everywhere function f(i::Int,L::Int)

	W = div(L,2)
	N = W*L 
	leads = ["A","B"]

	netT = Dict(l1*l2 => rand() for l1=leads for l2=leads if l1!=l2) 

	siteT = Dict(l=>rand(2,N) for l in leads)
	ldos  = rand(N) 
	dos = rand() 


	return Dict{String,Any}("i" => i,
													"siteT" => siteT,
													"ldos" => ldos,
													"netT" => netT,
													"dos" => dos,
													)


end 

@everywhere function init_data(L::Int)

	W = div(L,2)
	N = W*L 
	leads = ["A","B"]

	netT = Dict(l1*l2 => SharedArray{Float64}(1,NR_ENERGIES) for l1=leads for l2=leads if l1!=l2) 

	siteT = Dict(l=>SharedArray{Float64}(2,N,NR_ENERGIES) for l in leads)
	ldos  = SharedArray{Float64}(N,NR_ENERGIES)
	dos = SharedArray{Float64}(1,NR_ENERGIES) 

	i = SharedArray{Float64}(1,NR_ENERGIES)

	return Dict{String,Any}("i" => i,
													"siteT" => siteT,
													"ldos" => ldos,
													"netT" => netT,
													"dos" => dos,
													)

end 
@everywhere function fill_data!(data,i::Int,L::Int)

	W = div(L,2)
	N = W*L 
	leads = ["A","B"]

	setindex!(data["i"],i,i)

	
	for l1=leads,l2=leads 
		l1==l2 && continue 

		setindex!(data["netT"][l1*l2],rand(),i)

	end 

	for l in leads 

		copy!(selectdim(data["siteT"][l],1,i),rand(2,N)) 

	end 

	copy!(selectdim(data["ldos"],1,i),rand(N))

	setindex!(data["dos"],rand(),i)

	return 

end 





@everywhere sMB(x) = Base.summarysize(x)/2^20 

@everywhere function get_data(L::Int)
	
		Utils.RecursiveMerge(1:NR_ENERGIES;dim=2,parallel=nworkers()>1) do i 
		
			f(i,L)
		
		end;
		return 
	end 
		
	
function fill_data_all!(shdata,L::Int)
	@sync @distributed for i=1:NR_ENERGIES

		fill_data!(shdata, i, L)

	end 

end 



for L  = [2,60,100,3][1:1]

	println() 

	@testset "L=$L" begin 

#		@showtime get_data(L)

		shdata = @showtime init_data(L) 


		shdata2 = MeshInPlace.init_storage(f(1,L), NR_ENERGIES; parallel=true,shared=true)



		for (k,s) in shdata

			if s isa AbstractArray 
				@test s≈shdata2[k]

			elseif s isa AbstractDict 
			
				for (q,sq) in s 

					@test sq≈shdata2[k][q]

				end 

			end 

		end 

#		@showtime fill_data_all!(shdata, L)


		GC.gc() 

		@test true 
	
	end 

end 

nothing 

