import myLibs:ObservablesFromGF ,Algebra, GreensFunctions, QuantumMechanics

leads = ["A", "B"]

function LeadLayerSlicer(L::AbstractString, uc...)

	(L=="A" ? "LeftLead" : "RightLead",1), (Colon(),)

end  


widths = Dict(L=>rand(10:20) for L in leads)

VirtLeads = Dict( Symbol(LeadLayerSlicer(L)[1][1])=>Dict(:intercell=>[rand(ComplexF64,w,w)]) for (L,w) in widths)



#SE = Dict(L=>rand(ComplexF64,widths[L],widths[L]) for L in leads)




GF = Dict((L1,L2) => rand(ComplexF64,widths[L1],widths[L2]) for L1 in leads for L2 in leads)


CC = ObservablesFromGF.CaroliConductance


function G(args...; kwargs...)::Matrix

	Ls = filter(a-> a isa AbstractString, args)

	length(Ls)==1 && return GF[(Ls[1],Ls[1])]

	length(Ls)==2 && return GF[Tuple(Ls)]

	error() 

end 

function GA(args...; kwargs...)::Matrix

	Ls = filter(a-> a isa AbstractString, args)

	length(Ls)==1 && return GF[(Ls[1],Ls[1])]'

	length(Ls)==2 && return GF[(Ls[2],Ls[1])]'

	error() 

end 

SE = Dict(L=>GreensFunctions.SelfEn_fromGDecim(VirtLeads, LeadLayerSlicer, L, 1)(G) for L in leads)

function get_SE(L::AbstractString, args...)::Matrix

	SE[L]

end 


function get_SEsd(L::AbstractString)::Function 

	f3252(gr::Function)::Matrix = SE[L]

end 



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

@testset "One lead" begin

for L in leads 

	ccL = CC(G, get_SE, L)

	for nr_proj in 0:2 
	
		projectors = [identity for i=1:nr_proj]

		@test ccL≈CC(G, L, get_SEsd(L), 1, projectors...)
	
		@test ccL≈CC(G, get_SE, L, 1, projectors...)
	
		@test ccL≈CC(G, get_SE, L, projectors...)
	
		@test ccL≈CC(G, get_SE, L, 1, projectors...)

		@test ccL≈CC(G, SE[L], L, 1, projectors...)

	end 

end 

end  





@testset "Gr, two leads" begin

for L1 in leads 

	ccL1 = CC(G, get_SE, L1)

	for L2 in leads 

		ccL2 = CC(G, SE[L1], SE[L2], L1, L2, 1)

		@test ccL2≈CC(G, get_SE, L1, L2)

		@test ccL2≈CC(G, L1, L2, get_SEsd(L1), get_SEsd(L2), 1)
		
		@test ccL2≈CC(G, L1, L2, VirtLeads, LeadLayerSlicer, 1)
		
		L2==L1 && @test ccL1 ≈ ccL2 

	end 

end 

end  







@testset "Gr and Ga, two leads" begin

for L1 in leads 

	ccL1 = CC(G, get_SE, L1)

	for L2 in leads 

		ccGG = CC(G,GA, get_SE, L1, L2) 
		
		@test ccGG≈CC(G, get_SE, L1, L2) 

		@test ccGG≈CC(G,GA, SE[L1], SE[L2], L1, L2, 1) 

		@test ccGG≈CC(G, GA, L1, L2, get_SEsd(L1), get_SEsd(L2), 1)
		
		@test ccGG≈CC(G, GA, L1, L2, VirtLeads, LeadLayerSlicer, 1)

		L2==L1 && @test ccL1 ≈ ccGG

	end 

end 

end  





#===========================================================================#
#
# integrated current 
#
#---------------------------------------------------------------------------#

CCi = ObservablesFromGF.CaroliCurrent 


function fG(E::ComplexF64)::Function 

	d=  imag(E) 


	@assert abs(d)>1e-10 
	
	return d<0 ? G : GA

end 


FD = QuantumMechanics.FermiDirac 

FDH = QuantumMechanics.FermiDiracHole 



	


@testset "current" begin 

for L1 in leads, L2 in leads 
	
	i = CCi(fG, 0.01, FD(0.1), FDH(0.3)∘-, L1, L2, get_SEsd(L1), get_SEsd(L2), 1)

	for nr_d in 1:2 

		deltas = [0.01,-0.01][1:nr_d]
	
		@test i≈CCi(fG, deltas..., FD(0.1), FDH(0.3)∘-, L1, L2, get_SEsd(L1), get_SEsd(L2), 1)
		
		@test i≈CCi(fG, deltas..., FD(0.1), FDH(0.3)∘-, L1, L2, VirtLeads, LeadLayerSlicer, 1)
		
	end  

end 



end 




println() 



for L in leads
	
	@show L 

	U,i,dir = GreensFunctions.inter_i_dir(VirtLeads, LeadLayerSlicer, L)

	@show size(U) i dir 

	getse = GreensFunctions.SelfEn_fromGDecim(VirtLeads, LeadLayerSlicer, L)
	
	se = GreensFunctions.SelfEn_fromGDecim(G, VirtLeads, LeadLayerSlicer)(L)

	@show size(getse(G)) size(se)

	@assert se≈getse(G)

	println()

end 










