using myLibs:  BandStructure 



b = BandStructure.Diagonalize(x->rand(ComplexF64, 10,10), rand(3,2))


@show size(b["Energy"])
