import myLibs: ObservablesFromGF 

Gr(args...) = rand(ComplexF64,5,5)

Hoppings = [rand(ComplexF64,5,5) for i=1:10]

function SE_lead(args...) 
	e = rand(ComplexF64,5,5)
end 

Bonds = [(1,2),(2,3),(3,4),(4,5),(2,5)]
RBonds = [[rand(2),rand(2)] for b in Bonds]

bondT = ObservablesFromGF.BondTransmission(Gr, Gr,Hoppings, Bonds, SE_lead, "A", 1;GaGr_check="warn")
bondT = ObservablesFromGF.BondTransmission(Gr,Hoppings, Bonds, SE_lead, "A", 1;GaGr_check="warn")

@show size(bondT)


siteT = ObservablesFromGF.SiteTransmission(Gr, Gr, Hoppings, Bonds, RBonds, SE_lead, "A", 1;GaGr_check="warn",dim=1)
siteT = ObservablesFromGF.SiteTransmission(Gr, Hoppings, Bonds, RBonds, SE_lead, "A", 1;GaGr_check="warn",dim=1)

@show size(siteT)


