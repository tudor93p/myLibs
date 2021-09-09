using myLibs:  BandStructure, Utils, Operators

nr_at = 10
nr_orb = 3 
nr_wf = 17 
nr_k = 7 
size_H = nr_at*nr_orb 

dim = 2
dim2 = [2,1][dim]
atoms = rand([nr_at,2][[dim,dim2]]...)


function H(k...) 

	h = rand(ComplexF64, size_H, size_H) .- 0.5 

	return h+h'

end 


kwargs = (dim=dim,nr_at=nr_at,nr_orb=nr_orb,size_H=size_H)

operators = [[
							"LDOS",
							"X"
							], [
									Operators.LDOS(;kwargs...),
									Operators.Position(1, atoms; kwargs...)
									]]



kPoints, kTicks = Utils.PathConnect(rand([3,2][[dim,dim2]]...), nr_k; dim=dim)

@show size(kPoints) kTicks

#filename = nothing
filename = (x::String) -> joinpath("Data/Bands",x)

@show BandStructure.FoundFiles_Bands(filename, operators[1], "jld")
println()

b = BandStructure.Diagonalize(H, round.(kPoints,digits=3), filename; kTicks=kTicks, storemethod="jld", operators=operators, nr_bands=nr_wf, dim=2)
@time  b = BandStructure.Diagonalize(H, round.(kPoints,digits=3), filename; kTicks=kTicks, storemethod="jld", operators=operators, nr_bands=nr_wf, dim=2)

println()

@show BandStructure.FoundFiles_Bands(filename, operators[1], "jld")


println()
c = BandStructure.Read_Bands(filename, operators[1], "jld")
@time c = BandStructure.Read_Bands(filename, operators[1], "jld")
println()

@testset "calc+read" begin 

for (k,v) in pairs(b)

	@show k size(v) 

	@test isapprox(v,c[k])

	println()

end 

end 

println()


