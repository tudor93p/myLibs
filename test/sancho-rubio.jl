import myLibs: GreensFunctions, Lattices, H_Superconductor,TBmodel
import PyPlot, LinearAlgebra

PyPlot.close.(1:10)

latt = Lattices.SquareLattice()

Dir = 1

width = 4 

Lattices.Superlattice!(latt, setindex!(width*[1,1],1,Dir))

Lattices.KeepDim!(latt, Dir) 

@show Lattices.NrVecs(Lattices.PosAtoms(latt))
println.(eachcol(Lattices.PosAtoms(latt))) 


PyPlot.figure(1) 

atoms(n) = Lattices.Atoms_ManyUCs(latt; Ns=n)

for n in 0:3 

	PyPlot.scatter(eachrow(atoms(n))...,label=n)

end 

PyPlot.legend() 

hopp = H_Superconductor.SC_Domain((ChemicalPotential=0,),	[1]) 

inter = TBmodel.HoppingMatrix(atoms(0),atoms(1);hopp...)
intra = TBmodel.HoppingMatrix(atoms(0);hopp...)

gf(E) = GreensFunctions.GF_SanchoRubio(E+0.001im, intra, inter; target="+")["+"]


@show typeof(gf(rand()))
@show size(gf(rand())) 



function analytical_g(E::Float64,nu::Int)

	Vlead = 0 #chem pot 

	tx,ty = 1,1 # hopp 

	e = Vlead - 2ty*cos(pi*nu/(width+1))

	if abs(E-e) > 2tx 
		
		return (E-e)/(2tx^2)*(1 - sqrt(1 - 4tx^2/(E-e)^2)) + 0.001im

	else 


		return (E-e)/(2tx^2) - im*sqrt(1/tx^2 - (E-e)^2/(4tx^2))

	end


end 




ENERGIES = LinRange(-4,4,100)

Y = zeros(2,width,length(ENERGIES))

for (iE,E) in enumerate(ENERGIES)

	evals = LinearAlgebra.eigvals(gf(E)) 

	for (iF,F) in enumerate((real,imag))
		
		Y[iF,:,iE] = F(evals)

	end 

end 


for (iY,nr_fig) in enumerate([2,3])

	PyPlot.figure(nr_fig)

	for iy in 1:width 

		PyPlot.scatter(ENERGIES, Y[iY,iy,:], c="k")

		PyPlot.plot(ENERGIES, (real,imag)[iY](analytical_g.(ENERGIES,iy)), label=iy)
	end 

	PyPlot.legend()

end 


# plots as in
# J Comput Electron (2013) 12:203–231 
# The recursive Green’s function method for graphene
# Caio H. Lewenkopf · Eduardo R. Mucciolo






