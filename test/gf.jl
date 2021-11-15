import myLibs: GreensFunctions, LayeredSystem, Lattices, Algebra, Utils, H_Superconductor,TBmodel, Operators, ObservablesFromGF, BandStructure, ArrayOps
import PyPlot, JLD ,SparseArrays , LinearAlgebra

colors = hcat(["brown","red","coral","peru"],["gold","olive","forestgreen","lightseagreen"],["dodgerblue","midnightblue","darkviolet","deeppink"]) 




colors = [colors[i,j] for i in axes(colors,1) for j in [3:size(colors,2);1:2]]

#colors = reshape(transpose(colors),:)
#],1,axis=1).
#T.reshape(-1)


PyPlot.close.(1:10)

LENGTH = 51

P1 = (length = LENGTH, Barrier_height = 2.0, Barrier_width = 0.03, SCDW_phasediff = 0., SCDW_p = 2, SCDW_width = 0.005, SCDW_position = 0, SCpx_magnitude = 0.4, delta = 0.002, width = div(LENGTH,2), SCpy_magnitude = 0.4, Hopping=1.0, ChemicalPotential=2.0)

P2 = (Attached_Leads = "AB",)

P3 = [(ChemPot = 0.0, Hopping = 1.0, Label = "A", Direction = -1, Contact = (-1, 1), Lead_Width = div(P1.width,2), Lead_Coupling = 1.0, SCbasis = true), (ChemPot = 0.0, Hopping = 1.0, Label = "B", Direction = 1, Contact = (1, -1), Lead_Width = max(1,div(P1.width,2)-0), Lead_Coupling = 1.0, SCbasis = true)]

enlim = [-0.39493671,  0.39493671]

Dev = Lattices.SquareLattice()

NxNy = [P1[:length], P1[:width]]

Lattices.ShiftAtoms!(Dev, n=(1 .- NxNy)/2)


selectMainDim(v::AbstractVector)::Number = v[1] 

selectMainDim(m::AbstractMatrix)::AbstractVector = m[1,:] 


												


Lattices.Superlattice!(Dev, NxNy, Labels=selectMainDim)

Lattices.ReduceDim!(Dev)


DevAtoms = Lattices.PosAtoms(Dev) 

PyPlot.figure(4)

PyPlot.scatter(eachrow(DevAtoms)..., label="Device",alpha=0.3)



@show Lattices.NrVecs(DevAtoms)

SurfaceAtoms = Lattices.SurfaceAtoms(DevAtoms,1)

SurfaceAtoms2 = Lattices.filterAtoms_byNrBonds(-5:-1, DevAtoms, 1)


@assert isapprox(SurfaceAtoms,SurfaceAtoms2) 


@show Lattices.NrVecs(SurfaceAtoms)  

CornerAtoms = Lattices.filterAtoms_byNrBonds(2, DevAtoms, 1)


@show Lattices.NrVecs(CornerAtoms)  

qs = Utils.Quadrant(CornerAtoms, Algebra.Mean(CornerAtoms, 2), dim=2) 

#println.(eachcol(CornerAtoms))
#@show qs  

CornerAtoms = CornerAtoms[:,sortperm(qs)]

println.(eachcol(CornerAtoms))

#		2 ----------- 1
#		|      |      |     
#		|	---- : ----	|
#		|      |      |
#		3 ----------- 4


LeadLatts = map(enumerate(P3)) do (ilp,lead_params)
	println()
	(x,y) = lead_params[:Contact]


	t,(a,b) = if abs(x) == 1 
								
								(y,x>0 ? [4,1] : [3,2])

						elseif abs(y)==1 

								 (x,y>0 ? [2,1] : [3,4])

						end 

	start = (CornerAtoms[:,a]*(1-t) + CornerAtoms[:,b]*(1+t))/2

	@show start 


	d = lead_params[:Direction]::Int 

	sDir,Dir = sign(d),abs(d)

	L = Lattices.SquareLattice(string("Lead", sDir>0 ? "+" : "-", Dir))
	
	nxny = repeat([lead_params[:Lead_Width]], 2)

	nxny[Dir] = sDir 

	Lattices.Superlattice!(L, nxny)


	Lattices.KeepDim!(L, Dir) 



	@show Lattices.NrVecs(Lattices.PosAtoms(L))
	@show Lattices.PosAtoms(L)[:,1]


	shift = start - Algebra.Mean(Lattices.PosAtoms(L), 2) 

#	@show shift 
	
	Lattices.ShiftAtoms!(L, r=shift)

	@show Lattices.PosAtoms(L)[:,1] 



	lv = Lattices.LattVec(L)
	@assert length(lv)==2 
	@show lv[1] 

	Lattices.Align_toAtoms!(L, SurfaceAtoms)
	

	@show Lattices.PosAtoms(L)[:,1]  
	if ilp==1 
		PyPlot.figure(4)
		PyPlot.scatter(eachrow(SurfaceAtoms)...,label="Surf",alpha=0.3)

	end 

#	PyPlot.scatter(eachrow(Lattices.PosAtoms(L))..., label=lead_params[:Label])

	return L

end 


println()

inds_DevBonds = Algebra.get_Bonds(DevAtoms, 1; dim=2) 

Rs_DevBonds =  Algebra.bondRs_fromInds(inds_DevBonds, DevAtoms; dim=2)

@show size(inds_DevBonds) size(Rs_DevBonds)

leadlabels = [string(l) for l in P2[:Attached_Leads]]

@show leadlabels


println()

function zerogap(args...; kwargs...)::Matrix{ComplexF64}

	zeros(ComplexF64, 4, 4)

end 


coupling_Hparam = map(P3) do lp
	
	(Hopping = lp[:Lead_Coupling], SC_Gap = (zerogap, 0),	)

end 

lead_Hparam = map(P3) do lp 
	
	(Hopping = lp[:Hopping], 
	 ChemicalPotential = lp[:ChemPot], 
	 SC_Gap = (zerogap, 0),
								)

end 

function local_potential(ri, rj=ri)#dev_params::UODict)::Union{Float64,Function}

	w = P1.Barrier_width * P1.length 

	prefactor = P1.Barrier_height*pi*w 

	return prefactor * Algebra.Lorentzian(first(ri), w)

end

function eta(R)

	t = tanh( first(R) / (P1.SCDW_width * P1.length) )
	
	eta0 = [P1[:SCpx_magnitude], P1[:SCpy_magnitude]*1im]
	
	phase_factor = exp(0.5im*P1.SCDW_phasediff*pi*sign(t))
	
	return eta0 .* setindex!(ones(2), t, P1.SCDW_p) .* phase_factor 


end 

if false 
PyPlot.figure(3)

X= LinRange(-P1.length/2,P1.length/2,300)

PyPlot.plot(X, [real(eta(x)[1]) for x in X],label="Re x")
PyPlot.plot(X, [real(eta(x)[2]) for x in X],label="Re y")

PyPlot.plot(X, [imag(eta(x)[1]) for x in X],label="Im x") 
PyPlot.plot(X, [imag(eta(x)[2]) for x in X],label="Im y")

PyPlot.plot(X, local_potential.(X), label="V") 

PyPlot.ylim(-1.1,1.1)

PyPlot.legend() 
end 

function chiral_gap(ri::AbstractVector, rj::AbstractVector)::Matrix{ComplexF64}

	eta_x,eta_y = eta(ri/2 + rj/2)
	


	delta(u)::Bool = isapprox(u,0,atol=0.0001)
	

	psi = 0  

	(x, y)= ri-rj 
	
	D = zerogap()
	
	D[1,4]=-1im*eta_y*(-delta(y - 1) + delta(y + 1))*delta(x)/2 + (1im*eta_x*delta(x - 1)/2 - 1im*eta_x*delta(x + 1)/2 + psi*delta(x))*delta(y)
	
	D[2,3]=-1im*eta_y*(-delta(y - 1) + delta(y + 1))*delta(x)/2 + (1im*eta_x*delta(x - 1)/2 - 1im*eta_x*delta(x + 1)/2 - psi*delta(x))*delta(y)
	
	D[3,2]=-1im*(-delta(y - 1) + delta(y + 1))*delta(x)*conj(eta_y)/2 + (1im*delta(x - 1)*conj(eta_x)/2 - 1im*delta(x + 1)*conj(eta_x)/2 - delta(x)*conj(psi))*delta(y)
	
	D[4,1]=-1im*(-delta(y - 1) + delta(y + 1))*delta(x)*conj(eta_y)/2 + (1im*delta(x - 1)*conj(eta_x)/2 - 1im*delta(x + 1)*conj(eta_x)/2 + delta(x)*conj(psi))*delta(y)
	
	return D
		
end 
			


dev_Hparam = (Hopping = P1[:Hopping],
							ChemicalPotential = P1[:ChemicalPotential],
							LocalPotential = local_potential,#(dev_params),
							SC_Gap = (chiral_gap, 1))

dev_Hopping = H_Superconductor.SC_Domain(dev_Hparam,	[1])

@assert maximum(abs, dev_Hopping[:Hopping]([1,0],[0,0])) > 0  


@show maximum(abs, dev_Hparam.SC_Gap[1]([10,0],[11,0])) 


for i in 1:size(DevAtoms,2)-1

	at = (DevAtoms[:,i],DevAtoms[:,i+1])

	Delta = dev_Hparam.SC_Gap[1](at...)

	NZ = SparseArrays.findnz(SparseArrays.sparse(Delta))

	if false #all(!isempty, NZ) 
		println(at[1],"->",at[2],": ")
		
		println.(eachrow(Delta)) 
		println()
		
	end# @show Delta#NZ...)
    

end 




@show maximum(abs, lead_Hparam[1].SC_Gap[1]([10,0],[11,0])) 

DevH = TBmodel.Bloch_Hamilt(Lattices.NearbyUCs(Dev,1); dev_Hopping...)

fname(x::String) = joinpath("test/testData",x) 

@show size(DevH())

#println("First atoms: ",eachcol(DevAtoms[:,1:2])...)

#println(DevH()[1:4,1:4])
#println()

IJV2 = ([1, 5, 8, 9, 12, 2, 6, 7, 10, 11, 3, 6, 7, 10, 11, 4, 5, 8, 9, 12, 1, 4, 5, 13, 16, 2, 3, 6, 14, 15, 2, 3, 7, 14, 15, 1, 4, 8, 13, 16, 1, 4, 9, 13, 16, 17, 20, 2, 3, 10, 14, 15, 18, 19, 2, 3, 11, 14, 15, 18, 19, 1, 4, 12, 13, 16, 17, 20, 5, 8, 9, 12, 13, 21, 24, 6, 7, 10, 11, 14, 22, 23, 6, 7, 10, 11, 15, 22, 23, 5, 8, 9, 12, 16, 21, 24, 9, 12, 17, 21, 24, 25, 28, 10, 11, 18, 22, 23, 26, 27, 10, 11, 19, 22, 23, 26, 27, 9, 12, 20, 21, 24, 25, 28, 13, 16, 17, 20, 21, 29, 32, 14, 15, 18, 19, 22, 30, 31, 14, 15, 18, 19, 23, 30, 31, 13, 16, 17, 20, 24, 29, 32, 17, 20, 25, 29, 32, 18, 19, 26, 30, 31, 18, 19, 27, 30, 31, 17, 20, 28, 29, 32, 21, 24, 25, 28, 29, 22, 23, 26, 27, 30, 22, 23, 26, 27, 31, 21, 24, 25, 28, 32], [1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 10, 10, 10, 10, 10, 10, 10, 11, 11, 11, 11, 11, 11, 11, 12, 12, 12, 12, 12, 12, 12, 13, 13, 13, 13, 13, 13, 13, 14, 14, 14, 14, 14, 14, 14, 15, 15, 15, 15, 15, 15, 15, 16, 16, 16, 16, 16, 16, 16, 17, 17, 17, 17, 17, 17, 17, 18, 18, 18, 18, 18, 18, 18, 19, 19, 19, 19, 19, 19, 19, 20, 20, 20, 20, 20, 20, 20, 21, 21, 21, 21, 21, 21, 21, 22, 22, 22, 22, 22, 22, 22, 23, 23, 23, 23, 23, 23, 23, 24, 24, 24, 24, 24, 24, 24, 25, 25, 25, 25, 25, 26, 26, 26, 26, 26, 27, 27, 27, 27, 27, 28, 28, 28, 28, 28, 29, 29, 29, 29, 29, 30, 30, 30, 30, 30, 31, 31, 31, 31, 31, 32, 32, 32, 32, 32], ComplexF64[2.000006359300477 + 0.0im, 1.0 + 0.0im, -0.2 + 0.0im, 1.0 + 0.0im, 0.0 + 0.2im, 2.000006359300477 + 0.0im, 1.0 + 0.0im, -0.2 + 0.0im, 1.0 + 0.0im, 0.0 + 0.2im, -2.000006359300477 + 0.0im, 0.2 + 0.0im, -1.0 + 0.0im, 0.0 + 0.2im, -1.0 + 0.0im, -2.000006359300477 + 0.0im, 0.2 + 0.0im, -1.0 + 0.0im, 0.0 + 0.2im, -1.0 + 0.0im, 1.0 + 0.0im, 0.2 + 0.0im, 2.000006359300477 + 0.0im, 1.0 + 0.0im, 0.0 + 0.2im, 1.0 + 0.0im, 0.2 + 0.0im, 2.000006359300477 + 0.0im, 1.0 + 0.0im, 0.0 + 0.2im, -0.2 + 0.0im, -1.0 + 0.0im, -2.000006359300477 + 0.0im, 0.0 + 0.2im, -1.0 + 0.0im, -0.2 + 0.0im, -1.0 + 0.0im, -2.000006359300477 + 0.0im, 0.0 + 0.2im, -1.0 + 0.0im, 1.0 + 0.0im, 0.0 - 0.2im, 2.0000544629349473 + 0.0im, 1.0 + 0.0im, -0.2 + 0.0im, 1.0 + 0.0im, 0.0 + 0.2im, 1.0 + 0.0im, 0.0 - 0.2im, 2.0000544629349473 + 0.0im, 1.0 + 0.0im, -0.2 + 0.0im, 1.0 + 0.0im, 0.0 + 0.2im, 0.0 - 0.2im, -1.0 + 0.0im, -2.0000544629349473 + 0.0im, 0.2 + 0.0im, -1.0 + 0.0im, 0.0 + 0.2im, -1.0 + 0.0im, 0.0 - 0.2im, -1.0 + 0.0im, -2.0000544629349473 + 0.0im, 0.2 + 0.0im, -1.0 + 0.0im, 0.0 + 0.2im, -1.0 + 0.0im, 1.0 + 0.0im, 0.0 - 0.2im, 1.0 + 0.0im, 0.2 + 0.0im, 2.0000544629349473 + 0.0im, 1.0 + 0.0im, 0.0 + 0.2im, 1.0 + 0.0im, 0.0 - 0.2im, 1.0 + 0.0im, 0.2 + 0.0im, 2.0000544629349473 + 0.0im, 1.0 + 0.0im, 0.0 + 0.2im, 0.0 - 0.2im, -1.0 + 0.0im, -0.2 + 0.0im, -1.0 + 0.0im, -2.0000544629349473 + 0.0im, 0.0 + 0.2im, -1.0 + 0.0im, 0.0 - 0.2im, -1.0 + 0.0im, -0.2 + 0.0im, -1.0 + 0.0im, -2.0000544629349473 + 0.0im, 0.0 + 0.2im, -1.0 + 0.0im, 1.0 + 0.0im, 0.0 - 0.2im, 2.0000544629349473 + 0.0im, 1.0 + 0.0im, 0.2 + 0.0im, 1.0 + 0.0im, 0.0 + 0.2im, 1.0 + 0.0im, 0.0 - 0.2im, 2.0000544629349473 + 0.0im, 1.0 + 0.0im, 0.2 + 0.0im, 1.0 + 0.0im, 0.0 + 0.2im, 0.0 - 0.2im, -1.0 + 0.0im, -2.0000544629349473 + 0.0im, -0.2 + 0.0im, -1.0 + 0.0im, 0.0 + 0.2im, -1.0 + 0.0im, 0.0 - 0.2im, -1.0 + 0.0im, -2.0000544629349473 + 0.0im, -0.2 + 0.0im, -1.0 + 0.0im, 0.0 + 0.2im, -1.0 + 0.0im, 1.0 + 0.0im, 0.0 - 0.2im, 1.0 + 0.0im, -0.2 + 0.0im, 2.0000544629349473 + 0.0im, 1.0 + 0.0im, 0.0 + 0.2im, 1.0 + 0.0im, 0.0 - 0.2im, 1.0 + 0.0im, -0.2 + 0.0im, 2.0000544629349473 + 0.0im, 1.0 + 0.0im, 0.0 + 0.2im, 0.0 - 0.2im, -1.0 + 0.0im, 0.2 + 0.0im, -1.0 + 0.0im, -2.0000544629349473 + 0.0im, 0.0 + 0.2im, -1.0 + 0.0im, 0.0 - 0.2im, -1.0 + 0.0im, 0.2 + 0.0im, -1.0 + 0.0im, -2.0000544629349473 + 0.0im, 0.0 + 0.2im, -1.0 + 0.0im, 1.0 + 0.0im, 0.0 - 0.2im, 2.000006359300477 + 0.0im, 1.0 + 0.0im, 0.2 + 0.0im, 1.0 + 0.0im, 0.0 - 0.2im, 2.000006359300477 + 0.0im, 1.0 + 0.0im, 0.2 + 0.0im, 0.0 - 0.2im, -1.0 + 0.0im, -2.000006359300477 + 0.0im, -0.2 + 0.0im, -1.0 + 0.0im, 0.0 - 0.2im, -1.0 + 0.0im, -2.000006359300477 + 0.0im, -0.2 + 0.0im, -1.0 + 0.0im, 1.0 + 0.0im, 0.0 - 0.2im, 1.0 + 0.0im, -0.2 + 0.0im, 2.000006359300477 + 0.0im, 1.0 + 0.0im, 0.0 - 0.2im, 1.0 + 0.0im, -0.2 + 0.0im, 2.000006359300477 + 0.0im, 0.0 - 0.2im, -1.0 + 0.0im, 0.2 + 0.0im, -1.0 + 0.0im, -2.000006359300477 + 0.0im, 0.0 - 0.2im, -1.0 + 0.0im, 0.2 + 0.0im, -1.0 + 0.0im, -2.000006359300477 + 0.0im]);







if LENGTH==4 
	
		@show size.(IJV2)
		IJV1=  SparseArrays.findnz(DevH())
	
	for item in zip(IJV1,IJV2)
	
	@assert isapprox(item...)
	
	end 

end 



println(P1) 
println(P2)
println.(P3) 
println()
spectra = nothing 
#spectra = Dict()

if !isnothing(spectra)

	DevDiag = BandStructure.Diagonalize(DevH, 
																			Lattices.BrillouinZone(Dev),
																			fname,
																			dim=2) 
	
	
	
	E2 = sort(DevDiag["Energy"])[1:10]
	E1= [-4.629323740217191, -4.629323740217172, -3.638156713627021, -3.6381567136270054, -2.6333392035002894, -2.633339203500278, -2.4192701553998983, -2.4192701553998877, -1.6610804434322783, -1.6610804434322735]
	
	LENGTH==4&&@assert isapprox(E1,E2)
	
	
	@show extrema(DevDiag["kLabels"])
	
	println() 
	
	
	spectra["Device"] = (DevDiag["kLabels"],DevDiag["Energy"])

end 

leads = map(zip(leadlabels, LeadLatts, coupling_Hparam, lead_Hparam, P3)) do (
											label, latt, coupling_HP, lead_HP, lp)

	hopps = hopp_c, hopp_l = map([coupling_HP, lead_HP]) do Params 

		hopp = H_Superconductor.SC_Domain(Params,	[1])

		@assert hopp[:Hopping]([20.5, -9.5], [20.5, -8.5])[1]==Params.Hopping

		return  hopp 

	end 


	if !isnothing(spectra) && !haskey(spectra, "Lead")
		
		@show Lattices.BrillouinZone(latt)

		kPoints, = Utils.PathConnect(Lattices.BrillouinZone(latt), 100, dim=2)

		LeadH = TBmodel.Bloch_Hamilt(Lattices.NearbyUCs(latt,1); argH="k", hopp_l...)

		LeadDiag = BandStructure.Diagonalize(LeadH, kPoints, dim=2)

		spectra["Lead"] = (LeadDiag["kLabels"],LeadDiag["Energy"])

	end 
	

	hc, hl = map(hopps) do hopp

		function hm(args...) 
			
			TBmodel.HoppingMatrix(args...; hopp...)

		end 

	end 
	
	
	@assert hopp_l[:Nr_UCs]<=1 

  get_TBL(L::Lattices.Lattice) = Lattices.NearbyUCs(L, hopp_l[:Nr_UCs])

	intra, inter = TBmodel.BlochHamilt_ParallPerp(latt, get_TBL, hopp_l)

	# intra might be a function of k,
	# inter might be more than one matrix (m ∈ {1,2,...})

	function gf(E::Number; k=[])::Matrix{ComplexF64}

		GreensFunctions.GF_SanchoRubio(E, intra(k), only(values(inter));
																	 target="+")["+"]
		
	end

	return LayeredSystem.PrepareLead( label, latt, hc, hl, gf)

end 

if !isnothing(spectra) 

	PyPlot.figure(1) 
	
	PyPlot.scatter(spectra["Lead"]..., label="Lead")
	
	x,y = spectra["Device"]
	
	if isnothing(enlim)
	
		PyPlot.scatter(x,y,label="Device")
	
	else 
	
		inds = enlim[1] .<= y .<= enlim[2] 
	
		PyPlot.scatter(LinRange(0,1,count(inds)), y[inds], label="Device")
	
	
		PyPlot.ylim(enlim)
	
	end 
	
	PyPlot.xlim(0,1)
	
	PyPlot.legend() 

end 



PyPlot.figure(4) 

for lead in leads 
	PyPlot.scatter(eachrow(lead[:head][1])...,label=lead[:label])

end 

PyPlot.legend()
#PyPlot.close(4)


observables = ["QP-Caroli",
							 "QP-BondVectorTransmission",
								"QP-SiteVectorTransmission",]

@show keys(leads[1])

isBond = Algebra.EuclDistEquals(1; tol=1e-5,dim=2)




layers = parse.(Int, Lattices.labelToComponent(Dev)) .+ 1 

@show extrema(unique(layers))

	
la,ial = Utils.FindPartners(enumerate(layers), sortfirst=true)

default_LayerAtomRels = Dict(

			:NrLayers=> maximum(layers),

			:LayerOfAtom => la,

			:IndsAtomsOfLayer => ial,

			:AtomsOfLayer => function AtomsOfLayer(L::Int)::AbstractMatrix
			
				DevAtoms[:, ial(L)]

										end
							)



NewGeom = LayeredSystem.NewGeometry(DevAtoms, default_LayerAtomRels;
#																		"forced"; 
																		Leads=leads, isBond=isBond, dim=2,
																		dev_Hopping...)

@assert length(NewGeom)==4


println()

LayerAtom,Slicer,LeadRels,VirtLeads=NewGeom 


G = GreensFunctions.GF_Decimation(
					dev_Hopping, VirtLeads, Slicer;
					LayerAtom...,
					)





hoppings = LayeredSystem.get_hoppings(dev_Hopping[:Hopping], 
																			Slicer, 
																			VirtLeads)


BondHoppings = [hoppings(iR...)' for iR in zip(inds_DevBonds, Rs_DevBonds)]

#G, Energy, observables, Slicer, VirtLeads 

trace_bond = trace_atoms = LinearAlgebra.tr 

#trace_bond = Operators.Trace("atoms", ones(4); dim=2, nr_orb=4, sum_up=true, nr_at=1)
#trace_atoms = Operators.Trace("atoms", ones(4); dim=2, nr_orb=4, sum_up=true)




#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

Energy = (rand()-0.5)*0.35

@show Energy 

g_ = G(Energy + P1.delta*im) 

self_energy_ = GreensFunctions.SelfEn_fromGDecim(g_, VirtLeads, Slicer)







layer = rand(3:LayerAtom[:NrLayers]-3)

@show layer 

layer_right = layer+1 

layer_left = layer-1 

A0 = LayerAtom[:AtomsOfLayer](layer)   


inds_atoms_of_layer = LayerAtom[:IndsAtomsOfLayer](layer)


i_atom_in_layer = rand(axes(inds_atoms_of_layer,1))

i_atom = inds_atoms_of_layer[i_atom_in_layer]

R_atom = A0[:,i_atom_in_layer] 

@show i_atom_in_layer i_atom R_atom 




function print_transm(ij, rs, t::Real, out::AbstractDict=Dict())

	isapprox(t,0,atol=1e-14) && return 

	rs_ = map([1,2]) do k 
		
		rk = rs[k]

		if ij[k]==i_atom 
			
			@assert isapprox(rk,R_atom,atol=1e-5)

			return "R"#string(Int.(rk))

		end 

		dr = rk-rs[[2,1][k]]

		e = argmax(abs.(dr)) 

		return string("R",dr[e]>0 ? "+" : "-","xy"[e])

	end

	println("T",Tuple(ij),"=",round(t,digits=4),": ",join(rs_,"=>"))

	out[Tuple(ij)] = t 

end  




layer_siteT = zero(A0)



#	A = LayerAtom[:AtomsOfLayer](L)

#	H_lL = TBmodel.HoppingMatrix(A0, A; dev_Hopping...)

#	@show H_lL
	


Gi = g_("Layer", layer, "A", 3)


GWG = Gi*GreensFunctions.DecayWidth(self_energy_("A",3))*Gi'  



nr_orb = 4



out1=Dict()

nr_at = size(A0,2)

H_ll = TBmodel.HoppingMatrix(A0; dev_Hopping...)

println()

@info "New method"

println()

println("Intra-layer BondTij=GWG H")

for i in 1:nr_at, j in i+1:nr_at

	I = TBmodel.Hamilt_indices(1:nr_orb, i, nr_orb)  
	J = TBmodel.Hamilt_indices(1:nr_orb, j, nr_orb)  


	Tij = ObservablesFromGF.BondTij(GWG[I,J], H_ll[J,I])
																	
	ObservablesFromGF.addSiteTiTj!(layer_siteT, 2, 
																 (i,j), 
																 (A0[:,i], A0[:,j]),
																 Tij) 

	if i_atom_in_layer in [i,j]

		print_transm(inds_atoms_of_layer[[i,j]], (A0[:,i], A0[:,j]), Tij, out1)

	end 

end  






A_left = LayerAtom[:AtomsOfLayer](layer-1)
A_right = LayerAtom[:AtomsOfLayer](layer+1)



H_toleft = TBmodel.HoppingMatrix(A0,A_left; dev_Hopping...)
H_toright = TBmodel.HoppingMatrix(A0,A_right; dev_Hopping...)





g_left  = g_("Layer", layer-1, "Layer", layer-1)
g_right = g_("Layer", layer+1, "Layer", layer+1)


sigma_left = GreensFunctions.SelfEn(H_toleft', 
																		g_("Layer", layer-1; dir="left")) 

#sigma_left = GreensFunctions.SelfEn(H_toleft', g_left)

#sigma_right = GreensFunctions.SelfEn(H_toright',g_right)

sigma_right = GreensFunctions.SelfEn(H_toright',
																		 g_("Layer",layer+1; dir="right"))


CM5 = ObservablesFromGF.longitudinal_conductance_matrix(
					g_("Layer",layer),
					sigma_left,
					sigma_right)



#for (i,Ri) in enumerate(eachcol(A_left))

#	for (j,Rj) 


println("\nLayer $(layer-1)=>$layer: tr(conductance_matrix)")


for (i,j) in Algebra.get_Bonds(A_left, A0, isBond, dim=2)


	I = TBmodel.Hamilt_indices(1:nr_orb, i, nr_orb)

	J = TBmodel.Hamilt_indices(1:nr_orb, j, nr_orb)


	Tij = real(LinearAlgebra.tr(CM5[I,J])) 

	ObservablesFromGF.addSiteTiTj!(layer_siteT, 2, 
																 (j,), 
																 (A_left[:,i], A0[:,j]),
																 Tij)

	if j==i_atom_in_layer


		print_transm((LayerAtom[:IndsAtomsOfLayer](layer-1)[i],i_atom),
								 (A_left[:,i], A0[:,j]), Tij, out1)

	end 

#from layer-1 to layer 


end 


A_left = LayerAtom[:AtomsOfLayer](layer-1+1)
A_right = LayerAtom[:AtomsOfLayer](layer+1+1)


A0 = LayerAtom[:AtomsOfLayer](layer+1)

H_toleft = TBmodel.HoppingMatrix(A0, A_left; dev_Hopping...) 

H_toright = TBmodel.HoppingMatrix(A0, A_right; dev_Hopping...)





g_left  = g_("Layer", layer-1+1)
g_right = g_("Layer", layer+1+1)


CM6 = ObservablesFromGF.longitudinal_conductance_matrix(
					g_("Layer", layer+1, "Layer", layer+1),
					GreensFunctions.SelfEn(H_toleft', g_left),
					GreensFunctions.SelfEn(H_toright', g_right))





println("\nLayer $layer=>$(layer+1): tr(conductance_matrix)")


for (i,j) in Algebra.get_Bonds(A_left, A0, isBond, dim=2)

	I = TBmodel.Hamilt_indices(1:nr_orb, i, nr_orb)   

	J = TBmodel.Hamilt_indices(1:nr_orb, j, nr_orb)  

	Tij = real(LinearAlgebra.tr(CM6[I,J]))

	ObservablesFromGF.addSiteTiTj!(layer_siteT, 2, 
																 (i,), 
																 (A_left[:,i], A0[:,j]),
																 Tij)

	if i==i_atom_in_layer

		print_transm((i_atom, LayerAtom[:IndsAtomsOfLayer](layer+1)[j]),
								 (A_left[:,i], A0[:,j]), Tij, out1)
		
	end 

#from layer to layer+1


end 







println("\nSite T: ",round.(layer_siteT[:, i_atom_in_layer],digits=4),"\n")


























#gA(i)::AbstractMatrix{ComplexF64} = g_("Atom",i,"A",2)

bvt = ObservablesFromGF.BondTransmission(g_, 
																				 BondHoppings,
																	inds_DevBonds, 
																	self_energy_, "A",2;
																	f=trace_bond) 



mask = [i_atom in I for I in inds_DevBonds]

out2 = Dict() 

println()
@info "Old method"
println()

for (k,rs,v) in zip(inds_DevBonds[mask],Rs_DevBonds[mask], bvt[mask])

	print_transm(k, rs, v, out2)

	test = if haskey(out1,k)

		isapprox(out1[k],v,atol=1e-10)
	
				elseif haskey(out1,reverse(k))

		isapprox(out1[reverse(k)],-v,atol=1e-10)

				elseif isapprox(v,0,atol=1e-10) 

		true 

				else 

					error("Different bonds!\n",keys(out1))

	end 

	test && println("\t\t\tCoincide")

	!test && @warn "\t\tDifferent"

end 


PyPlot.close(4) 

println()

# G already contains BondHoppings VirtLeads Slicer; but not inds -- needed? 


#bvt2 = ObservablesFromGF.BondTransmission2(G,
#																		Energy + P1.delta*im,
#																		inds_DevBonds,
#																		"A";
#																		f=trace_bond)





#@show isapprox(bvt,bvt2,atol=1e-10)
#@assert isapprox(bvt,bvt2,atol=1e-10)

#error()



bvt0 = ObservablesFromGF.BondTransmission0(g_, 
																BondHoppings,
																	inds_DevBonds;
																	f=trace_bond) 



jx0_ = ObservablesFromGF.SiteTransmission0(g_, 
																BondHoppings,
																	inds_DevBonds, Rs_DevBonds;
																	f=trace_bond, dim=1)
#@show size(jx0_)

close_to_dw = isapprox.(DevAtoms[1,:], minimum(abs,DevAtoms[1,:]))  


jx = ObservablesFromGF.SiteTransmission_fromBondT(bvt, inds_DevBonds, Rs_DevBonds; dim=1)

println("Site transmission: ",round.(jx[i_atom,:],digits=4)) 

println()
PyPlot.close(4)
error()  

@show transpose(jx[LayerAtom[:IndsAtomsOfLayer](layer),:])

error() 


jx_ = ObservablesFromGF.SiteTransmission(g_, 
																				 BondHoppings,
																	inds_DevBonds, Rs_DevBonds,
																	self_energy_, "A",2;
																	f=trace_bond, dim=1)




@assert isapprox(jx,jx_)

PyPlot.figure(4)

#PyPlot.quiver(eachrow(DevAtoms)..., eachcol(jx)...; pivot="middle")

PyPlot.quiver(eachrow(hcat(DevAtoms,A0))..., eachrow(hcat(transpose(jx),siteT))...; pivot="middle")



#ObservablesFromGF.CaroliConductance2 

gSD = g_("A", 2, "B", 2)


sigmaS = GreensFunctions.SelfEn_fromGDecim(g_, VirtLeads, Slicer)("A", 2)
sigmaD = GreensFunctions.SelfEn_fromGDecim(g_, VirtLeads, Slicer)("B", 2)

cm = ObservablesFromGF.longitudinal_conductance_matrix(gSD, sigmaS, sigmaD)




println("separate Caroli:")
foreach(println∘real,
				[
				 ObservablesFromGF.CaroliConductance2(gSD, sigmaS, sigmaD),
				 ObservablesFromGF.CaroliConductance(gSD, sigmaS, sigmaD),
				 ObservablesFromGF.CaroliConductance(g_, GreensFunctions.SelfEn_fromGDecim(g_, VirtLeads, Slicer), "A", "B", 2),
				 LinearAlgebra.tr(cm),
				 ])
println()

lA = leads[1] 

@assert lA[:label]=="A"


atoms2 = lA[:head][1] 

#atoms1 = Lattices.Atoms_ManyUCs(LeadLatts[1],Ns=1)

@assert isapprox(atoms2, Lattices.Atoms_ManyUCs(LeadLatts[1],Ns=0))

atoms(uc::Int=0) = Lattices.Atoms_ManyUCs(LeadLatts[1]; Ns=uc)  



#tr_orb = Operators.Trace(:orbitals; dim=2, nr_orb=4, nr_at=size(atoms1,2))


#@assert all(<(1e-12),imag(tr_orb(cm))) 
#
#lead_contacts = LayeredSystem.get_LeadContacts(DevAtoms; Leads=[lA], isBond=isBond)[1]
#
#
#
#atoms0 = DevAtoms[:,lead_contacts]
#
#gSD2 = g_("A", 1, "B", 2)
#sigmaS2 = GreensFunctions.SelfEn_fromGDecim(g_, VirtLeads, Slicer)("A", 1) 
#sigmaD2 = GreensFunctions.SelfEn_fromGDecim(g_, VirtLeads, Slicer)("B", 2)
#
#cm2 = real(tr_orb(ObservablesFromGF.longitudinal_conductance_matrix(gSD2, sigmaS2, sigmaD2)))



lead_intra_bond_inds = Algebra.get_Bonds(atoms(), isBond) 

lead_inter_bond_inds = Algebra.get_Bonds(atoms(0), atoms(1), isBond) 

uc=3  # unit cell of interest 


#siteT = zeros(2, size(atoms(uc),2))





for i in 1:7 

	@assert isapprox(self_energy_("A",i),self_energy_("A",uc))

end 




neighb_uc = uc#in [uc-1,uc+1]


bondT = SparseArrays.spzeros(Float64, 
														 size(atoms(uc),2), 
														 size(atoms(neighb_uc),2))



Gi = g_("A", uc, "A", neighb_uc)

A = Gi*GreensFunctions.DecayWidth(self_energy_("A",neighb_uc))*Gi' 


B = lA[:intracell][1]  

#only i>j  should be considered 

nr_orb = 4

for j in axes(bondT,2)

	J = TBmodel.Hamilt_indices(1:nr_orb, j, nr_orb)  

	for i in axes(bondT,1)

		I = TBmodel.Hamilt_indices(1:nr_orb, i, nr_orb) 

		Tij = ObservablesFromGF.BondTij(view(A,I,J), view(B,J,I))

#		transversal formula gives nonzero for lead intra-cell 
#		transversal formula gives zeros for lead inter-cell 



		isapprox(Tij,0,atol=1e-14) || setindex!(bondT, Tij, i, j)

	end 

end 

@show SparseArrays.findnz(bondT)


neighb_uc = uc-1  # ==drain. Source==uc

bondTl = zeros(Float64, length(lead_inter_bond_inds))

CM = ObservablesFromGF.longitudinal_conductance_matrix(
															 g_("A",uc, "A", neighb_uc),
															 self_energy_("A",uc),
															 GreensFunctions.SelfEn(lA[:intercell][1],
																											g_("A",neighb_uc)))


for (bond_index,(i,j)) in enumerate(lead_inter_bond_inds)

	I = TBmodel.Hamilt_indices(1:nr_orb, i, nr_orb) 
	J = TBmodel.Hamilt_indices(1:nr_orb, j, nr_orb)  

	Tij = LinearAlgebra.tr(CM[I,J])

	@assert isapprox(imag(Tij),0,atol=1e-14)

	bondTl[bond_index] = real(Tij)

end 



#uc=>uc: bondT 
#uc=>uc-1: bondTl 
#uc-1=>uc: bondTl 



siteT = zeros(Float64,size(atoms(uc)))


A0 = atoms(uc)
AL = atoms(uc-1)
AR = atoms(uc+1) 

nr_at = size(A0,2)


for i in 1:nr_at, j in i+1:nr_at 

	ObservablesFromGF.addSiteTiTj!(siteT, 2, (i,j), (A0[:,i],A0[:,j]), bondT[i,j])

end 

for (bond_index,(i,j)) in enumerate(lead_inter_bond_inds) 

	ObservablesFromGF.addSiteTiTj!(siteT, 2, (i,j), 
																 (A0[:,i], AL[:,j]), bondTl[bond_index])


end 



#PyPlot.quiver(eachrow(A0)..., eachrow(siteT)...; pivot="middle")










error()



#jA = ObservablesFromGF.SiteTransmission(g_, 
#																				 BondHoppings,
#																	inds_DevBonds, Rs_DevBonds,
#																	self_energy_, "A",2;
#																	f=trace_bond, dim=1)
#
#
##ObservablesFromGF.SiteTransmission_fromBondT(real(tr_orb(cm)),
#
#
#PyPlot.figure(4) 
#
#PyPlot.quiver(eachrow(atoms1/2+atoms2/2)..., 
#							eachrow(ArrayOps.multiply_elementwise(real(tr_orb(cm)),atoms2-atoms1))...;
#							pivot="middle")
#
#
#PyPlot.quiver(eachrow(atoms1/2+atoms0/2)..., 
#							eachrow(ArrayOps.multiply_elementwise(cm2,atoms1-atoms0))...;
#							pivot="middle")
#
##PyPlot.close.(1:10)
error()





jx = jx[close_to_dw,1] 


jx0 = ObservablesFromGF.SiteTransmission_fromBondT(bvt0, inds_DevBonds, Rs_DevBonds; dim=1)



@assert isapprox(jx0,jx0_) 

jx0 = jx0[close_to_dw,1] 

@show unique(DevAtoms[1,close_to_dw])

@show size(jx) 



#@show ObservablesFromGF.TunnelingConductance_LeeFisher(g_, "A", 2)
#@show ObservablesFromGF.TunnelingConductance_LeeFisher(g_, "B", 2)








PyPlot.figure(5)
y = DevAtoms[2,close_to_dw]

PyPlot.plot(y,jx, label="j", c=colors[1])
PyPlot.plot(y,jx0, label="j0", c=colors[2])

PyPlot.xlabel("y")
PyPlot.ylabel("jx")

PyPlot.legend()









ENERGIES = LinRange(0,enlim[2],40)|>collect 

iE0 = argmin(abs.(ENERGIES.-abs(Energy)))



out = Utils.RecursiveMerge(collect(enumerate(ENERGIES)); dim=2) do (iE,En)

#	print("\rEnergy=",round(Energy,digits=3),"      ")

	g = G(En + P1.delta*im)


	
	d = Dict()
#	@show size(g("Atom",1,"Atom",2))

	self_energy = GreensFunctions.SelfEn_fromGDecim(g, VirtLeads, Slicer)

	l,k = leadlabels

#	println(maximum(abs, g(l,2,k,2) ))

	d["Energy"] = En

	uc = 5

	d["QP-Caroli"] = ObservablesFromGF.CaroliConductance(
												g(l, uc, k, uc), 
												self_energy(l, uc), 
												self_energy(k, uc);
												f=trace_atoms,
												)

	@assert isapprox(d["QP-Caroli"],
									 ObservablesFromGF.CaroliConductance(g, self_energy,
																											 l, k, uc))


	uc = 2

	d["QP-Caroli2"] = ObservablesFromGF.CaroliConductance2(g, self_energy,
																												 l, k, uc) 

	if iE==iE0 

		@show iE Energy  
	
		for (k,v) in pairs(d)
			println(k,": ",real(v))
		end 

	end 
#	println()
#@time ObservablesFromGF.CaroliConductance(g, self_energy, l, k, 2)
#	@time ObservablesFromGF.CaroliConductance2(g, self_energy, l, k, 2)

	return d 

end  

#println("\r                             ")


#@assert isapprox(out["QP-Caroli"], out["QP-Caroli2"]) 
println()


@show keys(out) 

@show size.(values(out))

for (k,v) in out 

	println(k,extrema(real(v)))

end


#@show [maximum(abs, v) for v in values(out)] 

println() 


xlabel = "QP-Caroli"



PyPlot.figure(2)


for (i,x) in enumerate(["QP-Caroli","QP-Caroli2"])
	
	for s in [-1,+1]

		@assert maximum(abs, imag(out[x]))<1e-12
	
		label = s==1 ? x : nothing  
	
		PyPlot.plot(real(out[x]),s*real(out["Energy"]),c=colors[i],label=label)
	end 

end 

PyPlot.plot([0,1.5],[Energy,Energy], alpha=0.4)

PyPlot.xlabel(xlabel);#PyPlot.ylabel(ylabel)


isnothing(enlim) || PyPlot.ylim(enlim)

PyPlot.xlim(0,1.5)

PyPlot.legend()


















