import myLibs: GreensFunctions, LayeredSystem, Lattices, Algebra, Utils, H_Superconductor,TBmodel, Operators, ObservablesFromGF, BandStructure
import PyPlot, JLD ,SparseArrays , LinearAlgebra

colors = hcat(["brown","red","coral","peru"],["gold","olive","forestgreen","lightseagreen"],["dodgerblue","midnightblue","darkviolet","deeppink"]) 




colors = [colors[i,j] for i in axes(colors,1) for j in [3:size(colors,2);1:2]]

#colors = reshape(transpose(colors),:)
#],1,axis=1).
#T.reshape(-1)


PyPlot.close.(1:10)

LENGTH = 60

P1 = (length = LENGTH, Barrier_height = 2.0, Barrier_width = 0.03, SCDW_phasediff = 0., SCDW_p = 2, SCDW_width = 0.005, SCDW_position = 0, SCpx_magnitude = 0.4, delta = 0.002, width = div(LENGTH,2), SCpy_magnitude = 0.4, Hopping=1.0, ChemicalPotential=2.0)

P2 = (Attached_Leads = "AB",)

P3 = [(ChemPot = 0.0, Hopping = 1.0, Label = "A", Direction = -1, Contact = (-1, 1), Lead_Width = div(P1.width,2), Lead_Coupling = 1.0, SCbasis = true), (ChemPot = 0.0, Hopping = 1.0, Label = "B", Direction = 1, Contact = (1, -1), Lead_Width = div(P1.width,2), Lead_Coupling = 1.0, SCbasis = true)]

enlim = [-0.39493671,  0.39493671]

Dev = Lattices.SquareLattice()

NxNy = [P1[:length], P1[:width]]

Lattices.ShiftAtoms!(Dev, n=(1 .- NxNy)/2)


selectMainDim(v::AbstractVector)::Number = v[1] 

selectMainDim(m::AbstractMatrix)::AbstractVector = m[1,:] 


												


Lattices.Superlattice!(Dev, NxNy, Labels=selectMainDim)

Lattices.ReduceDim!(Dev)


DevAtoms = Lattices.PosAtoms(Dev) 

#PyPlot.figure(4)
#
#PyPlot.scatter(eachrow(DevAtoms)..., label="Device")



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


LeadLatts = map(P3) do lead_params 
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

#	PyPlot.figure(4)
#	PyPlot.scatter(eachrow(SurfaceAtoms)...,label="Surf")

	lv = Lattices.LattVec(L)
	@assert length(lv)==2 
	@show lv[1] 

	Lattices.Align_toAtoms!(L, SurfaceAtoms)
	

	@show Lattices.PosAtoms(L)[:,1]  #Lattices.PosAtoms(L)[:,2]

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

fname(x::String) = joinpath("testData",x) 

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
PyPlot.close(4)


observables = ["QP-Caroli",
							 "QP-BondVectorTransmission",
								"QP-SiteVectorTransmission",]

@show keys(leads[1])

isBond = Algebra.EuclDistEquals(1; tol=1e-5,dim=2)




layers = parse.(Int, Lattices.labelToComponent(Dev)) .+ 1 

@show unique(layers)

	
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





hoppings = LayeredSystem.get_hoppings(dev_Hopping[:Hopping], Slicer, VirtLeads)


BondHoppings = [hoppings(iR...)' for iR in zip(inds_DevBonds, Rs_DevBonds)]

#G, Energy, observables, Slicer, VirtLeads 

trace_bond = trace_atoms = LinearAlgebra.tr 

#trace_bond = Operators.Trace("atoms", ones(4); dim=2, nr_orb=4, sum_up=true, nr_at=1)
#trace_atoms = Operators.Trace("atoms", ones(4); dim=2, nr_orb=4, sum_up=true)

Energy = (rand()-0.5)*0.3

@show Energy 

g_ = G(Energy + P1.delta*im) 

self_energy_ = GreensFunctions.SelfEn_fromGDecim(g_, VirtLeads, Slicer)


#gA(i)::AbstractMatrix{ComplexF64} = g_("Atom",i,"A",2)

bvt = ObservablesFromGF.BondTransmission(g_, 
																				 BondHoppings,
																	inds_DevBonds, 
																	self_energy_, "A",2;
																	f=trace_bond) 

bvt0 = ObservablesFromGF.BondTransmission0(g_, 
																BondHoppings,
																	inds_DevBonds;
																	f=trace_bond) 


jx0_ = ObservablesFromGF.SiteTransmission0(g_, 
																BondHoppings,
																	inds_DevBonds, Rs_DevBonds;
																	f=trace_bond, dim=1)


close_to_dw = isapprox.(DevAtoms[1,:], minimum(abs,DevAtoms[1,:]))  


jx = ObservablesFromGF.SiteTransmission_fromBondT(bvt, inds_DevBonds, Rs_DevBonds; dim=1)

jx_ = ObservablesFromGF.SiteTransmission(g_, 
																				 BondHoppings,
																	inds_DevBonds, Rs_DevBonds,
																	self_energy_, "A",2;
																	f=trace_bond, dim=1)
@assert isapprox(jx,jx_)

jx = jx[close_to_dw,1] 


jx0 = ObservablesFromGF.SiteTransmission_fromBondT(bvt0, inds_DevBonds, Rs_DevBonds; dim=1)



@assert isapprox(jx0,jx0_) 

jx0 = jx0[close_to_dw,1] 

@show unique(DevAtoms[1,close_to_dw])

@show size(jx) 



@show ObservablesFromGF.TunnelingConductance_LeeFisher(g_, "A", 2)
@show ObservablesFromGF.TunnelingConductance_LeeFisher(g_, "B", 2)








PyPlot.figure(5)
y = DevAtoms[2,close_to_dw]

PyPlot.plot(y,jx, label="j", c=colors[1])
PyPlot.plot(y,jx0, label="j0", c=colors[2])

PyPlot.xlabel("y")
PyPlot.ylabel("jx")

PyPlot.legend()










ENERGIES = LinRange(0,enlim[2],30)|>collect 

out = Utils.RecursiveMerge(ENERGIES; dim=2) do Energy 

	print("\rEnergy=",round(Energy,digits=3),"      ")

	g = G(Energy + P1.delta*im)

	
	d = Dict()
#	@show size(g("Atom",1,"Atom",2))

	self_energy = GreensFunctions.SelfEn_fromGDecim(g, VirtLeads, Slicer)

	k,l = leadlabels

#	println(maximum(abs, g(l,2,k,2) ))

	d["Energy"] = Energy

	d["QP-Caroli"] = ObservablesFromGF.CaroliConductance(
												g(l, 2, k, 2), 
												self_energy(l, 2), 
												self_energy(k, 2);
												f=trace_atoms,
												)

	@assert isapprox(d["QP-Caroli"],
									 ObservablesFromGF.CaroliConductance(g, self_energy,
																											 l, k, 2))

	return d 

end  

println("\r                             ")



@show keys(out) 

@show size.(values(out))

for (k,v) in out 

	println(k,extrema(real(v)))

end


#@show [maximum(abs, v) for v in values(out)] 

println() 


xlabel = "QP-Caroli"



PyPlot.figure(2)


for (i,x) in enumerate(["QP-Caroli"])
	
	for s in [-1,+1]

		@assert maximum(abs, imag(out[x]))<1e-12
	
		label = s==1 ? x : nothing  
	
		PyPlot.plot(real(out[x]),s*real(out["Energy"]),c=colors[i],label=label)
	end 

end 

PyPlot.xlabel(xlabel);PyPlot.ylabel(ylabel)


isnothing(enlim) || PyPlot.ylim(enlim)

PyPlot.xlim(0,1.5)

PyPlot.legend()



















