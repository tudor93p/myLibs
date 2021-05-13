using myLibs:Lattices, LayeredSystem 


atoms = Lattices.PosAtoms(Lattices.Superlattice(Lattices.SquareLattice(), [10,7]))

isBond = Algebra.EuclDistEquals(1.0; dim=2) 


D = LayeredSystem.LayerAtomRels(atoms, "forced"; isBond=isBond, dim=2) 


@show D[:NrLayers]
@show D[:LayerOfAtom](10)
@show D[:IndsAtomsOfLayer](3)
@show D[:AtomsOfLayer](3)







