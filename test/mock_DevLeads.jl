using myLibs: Lattices, GreensFunctions, ArrayOps

import LinearAlgebra

function get_hopp(nr_orb_)

	intercell = rand(ComplexF64,nr_orb_,nr_orb_)  # h(ri,rj) != h(rj,ri)'

	intercell .+= intercell'


	function hopp(ri,rj)

		local_pot = ArrayOps.UnitMatrix(nr_orb_)*isapprox(LinearAlgebra.norm(ri-rj),0,atol=1e-8)*2.0  
		
		
		d1 = isapprox(LinearAlgebra.norm(ri-rj), 1, atol=1e-8) 
	
		atom_hopp = intercell*d1  

		A = ArrayOps.UnitMatrix(nr_orb_, ComplexF64)

		for i=1:nr_orb_ 
			
			A[i,i] *= (-1)^i * im * (sum(ri)-sum(rj))

		end 

		A.*=d1 

		return local_pot + atom_hopp + A 
	
	end  

end 

function get_latt(Nxy::Vector{Int})

	l1 = Lattices.SquareLattice()

	Lattices.ShiftAtoms!(l1, n=(1 .- Nxy)/2)

	Lattices.Superlattice!(l1, Nxy, Labels=first)  

	Lattices.ReduceDim!(l1)

	return l1 

end 

function get_atoms(Nxy::Vector{Int})

	Lattices.PosAtoms(get_latt(Nxy))

end 

get_leadlatt(n::Int=1, label...) = Lattices.KeepDim!(
																					Lattices.Superlattice( 
																					 Lattices.SquareLattice(label...),
																					 [1,n]
																					 )
																					 , 1)
get_lead(L::Lattices.Lattice, args...) = GreensFunctions.PrepareLead(only(Lattices.sublatt_labels(L)), L, args...)

get_lead(n::Int, label::String, args...) = GreensFunctions.PrepareLead(label,get_leadlatt(n, label), args...)


d_nn = LinearAlgebra.norm(diff(get_atoms([2,1]),dims=2))


isBond = Algebra.EuclDistEquals(d_nn; dim=2)
 



function default_LayerAtomRels(DL::Lattices.Lattice)

	atoms = Lattices.PosAtoms(DL)
	layers = parse.(Int, Lattices.labelToComponent(DL)) .+ 1

		# layer of each atom [1, 2, ..., nr_at] => [L(1), L(2), ..., L(nr_at)]
	
	la,ial = Utils.FindPartners(enumerate(layers), sortfirst=true)

	return Dict(

			:NrLayers=> maximum(layers),

			:LayerOfAtom => la,

			:IndsAtomsOfLayer => ial,

			:AtomsOfLayer => function AtomsOfLayer(L::Int)::AbstractMatrix
			
									Lattices.Vecs(atoms, ial(L))

										end
							)
end



function get_two_leads(
											 lead_width::Int,
											 atoms::AbstractMatrix,
											 labels::AbstractVector{String}=["A","B"];
											 hopp...
											 )
		
	[get_lead(Lattices.Align_toAtoms!(get_leadlatt(lead_width,lab),atoms,d),hopp,0.001) for (d,lab)=zip([-1,1],labels)]


end  


@testset "valid hopping" begin 

	Rs = rand(2,10) 
	for n=1:4 

		h = get_hopp(n)
	
		for ri=eachcol(Rs) 
	
			@test LinearAlgebra.ishermitian(h(ri,ri)) 

			for rj=eachcol(Rs)

				@test h(ri,rj)≈h(rj,ri)' 

				LinearAlgebra.norm(ri-rj)>1.1 && @test LinearAlgebra.norm(h(ri,rj))<1e-10

			end 
	
		end 

	end 

end 
