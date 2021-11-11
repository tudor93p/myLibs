###############################################################################
# This file refers to the model in the seminal paper of Kane and Mele
#   PRL 95, 146802 (2005)  
#   "Z2 Topological Order and the Quantum Spin Hall Effect"
###############################################################################

module H_Graphene



import ..LA 

import ..Utils, ..TBmodel, ..Algebra, ..BandStructure


export Graphene_Sheet#,Dirac_Points


println("Loaded H_Graphene")  



#===========================================================================#
#
# Describes hoppings between two atoms in different sheets
#
#---------------------------------------------------------------------------#

function Coupling_Sheets(param_H_,pyLatt,dist_tol=1e-5,hopp_cutoff=1e-6)#,partialH=false)

  nr_uc = 1 		# minimum, provided the kind of hoppings
  nr_neighbors = 2
  slcoord = Dict("A"=> 1.0, "B"=> -1.0)
  d0 = 1

  aux, TB_Input = pyLatt.TB_Output(nr_uc,nr_neighbors,dist_tol,slcoord)


  ms, UCs, AtomsUC = TB_Input 

        # ---------------- initialize parameters --------------------- #      


  default_param_H = (
	Interlayer_Hopping	= 0.0,
	Interlayer_HoppDecay	= 1.0/50.0,
	Coupled_Zs		= rand(2),
					)


  issubset(keys(param_H_),keys(default_param_H)) || error("The parameters must be among "*join(map(string,keys(default_param_H)),", ")*".")

  param_H = merge(default_param_H,param_H_)

  same = Utils.fSame(dist_tol)



        # ------------------- needed of the hoppings ----------------- #      


  Append,Sum,Hoppings = TBmodel.Add_Hopping_Terms(d0,hopp_cutoff)
		# value, condition, (function)

  maxdist,update = TBmodel.Update_Maxdist(dist_tol,hopp_cutoff)


        # ------------------- interlayer hopping --------------------- #     
 
  Zsep = param_H.Coupled_Zs |> x->maximum(x)-minimum(x)

  function interlayer(ri,rj)

    dist_xyz = LA.norm(ri[1:3].-rj[1:3]) 

    return abs2(ri[3]-rj[3])/dist_xyz^2*
                   exp(-(dist_xyz-Zsep)/param_H.Interlayer_HoppDecay)
  end


  Hoppings = Append(Hoppings, param_H[:Interlayer_Hopping],
		true,
		interlayer,
							)

        # -------------- compute max_dist and sum up ------------- #      

  function maxdist_() 

    Zs = param_H.Coupled_Zs; r0 = [0,0,Zs[1]]
  
    for d in range(0,1000Zsep,step=Zsep/20)
  
      isapprox(interlayer(r0,[d,0,Zs[2]]),0.0,atol=hopp_cutoff) && return d

    end
  end

  maxdist = update(maxdist,param_H[:Interlayer_Hopping],maxdist_())


  couple(ri,rj) = (same([ri[3],rj[3]],param_H.Coupled_Zs) || same([rj[3],ri[3]],param_H.Coupled_Zs)) &&  LA.norm(ri[1:2].-rj[1:2]) < maxdist




  At1 = AtomsUC[abs.(AtomsUC[:,3].-param_H.Coupled_Zs[1]).<dist_tol,1:2]

  if size(At1,1) > 0

    nr = count( Algebra.OuterDist(At1 .- At1[1:1,:], -UCs[:,1:2]) .< maxdist)
  else
  
   nr = 0

  end


  return TB_Input, (Sum(Hoppings,couple), nr, hopp_cutoff)



   
end



#===========================================================================#
#
# Describes hoppings between two atoms in the same sheet
#
#---------------------------------------------------------------------------#

function Sheet(param_H_,pyLatt,dist_tol=1e-5,hopp_cutoff=1e-6)#,partialH=false)

  nr_uc = 1 		# minimum, provided the kind of hoppings
  nr_neighbors = 2
  slcoord = Dict("A"=> 1.0, "B"=> -1.0)
  d0 = 1

  (dim, dist), TB_Input = pyLatt.TB_Output(nr_uc,nr_neighbors,dist_tol,slcoord)

  ms, UCs, AtomsUC = TB_Input # it happens that they're also needed here

  sl = size(AtomsUC,2)	# position where the sublattice is stored

  AllAtoms = Algebra.FlatOuterSum(AtomsUC,UCs)






        # ---------------- initialize parameters --------------------- #      

  default_param_H = (
	Intralayer_Hopping	= 1.0,
	Electric_Field_x	= 0.0,
	Electric_Field_y	= 0.0,
	Electric_Field_z	= 0.0,
	Lattice_Imbalance	= 0.0,
#        KaneMele_SOC		= 0.0,
#	Rashba_SOC		= 0.0,
	AntiHaldane		= 0.0,
					)


  param_H = Utils.Combine_NamedTuples(param_H_,default_param_H)



  same = Utils.fSame(dist_tol)

        # ------------------- needed of the hoppings ----------------- #      


  Append,Sum,Hoppings = TBmodel.Add_Hopping_Terms(d0,hopp_cutoff)
		# value, condition, (function)
  
  maxdist,update = TBmodel.Update_Maxdist(dist_tol,hopp_cutoff)

  
        # -------------- nearest neighbor ------------------------ #      

  Hoppings = Append(Hoppings, param_H[:Intralayer_Hopping],
#		(ri,rj) -> same(ri[1:2].-rj[1:2],dist[1]), REWRITE 
							)
  maxdist = update(maxdist,param_H[:Intralayer_Hopping],dist[1])
  
#  println();println(param_H[:Intralayer_Hopping]);println();

        # -------------- atom bias (lattice imbalance) ----------- #      

  Hoppings = Append(Hoppings, param_H[:Lattice_Imbalance],
#		(ri,rj) -> same(ri,rj), #REWRITE
#		(ri,rj) -> ri[sl], #REWRITE
							)

        # ------------------- electric field -------------------- #      

  Hoppings = Append(Hoppings, param_H[:Electric_Field_x],
#		(ri,rj) -> same(ri,rj), #REWRITE
#		(ri,rj) -> ri[1], #REWRITE
							)


  Hoppings = Append(Hoppings, param_H[:Electric_Field_y],
#		(ri,rj) -> same(ri,rj),
#		(ri,rj) -> ri[2],
							)


  Hoppings = Append(Hoppings, param_H[:Electric_Field_z]*(dim==3),
#		(ri,rj) -> same(ri,rj),
#		(ri,rj) -> ri[3],
							)

        # -------------- anti-Haldane nnn imaginary hopping ------ #      


  nn(ri) = findall(same.(Algebra.FlatOuterDist(hcat(ri[1:2]...),AllAtoms[:,1:2]),dist[1]))

  chirality(ri,rj) = sign(Algebra.cross(eachcol(hcat(ri,rj) .- AllAtoms[intersect(nn(ri),nn(rj))[1],:] )...)[3])


  Hoppings = Append(Hoppings, param_H[:AntiHaldane],
		(ri,rj) -> same(ri[1:2].-rj[1:2],dist[2]),
		(ri,rj) -> 1im/sqrt(27)*ri[sl]*chirality(ri,rj),
							)

  maxdist = update(maxdist,param_H[:AntiHaldane],dist[2])

        # -------------- Kane-Mele spin-orbit coupling ----------- #      

#    if same(la.norm(dif[:2]),dist_nnn)*same(abs(dif[2]),0.)*neighbors.size == 0
#  push!(hoppings, TBmodel.Hopping_Term( 
#		(ri,rj) -> isapprox(LA.norm(ri.-rj),dist_nnn,atol=dist_tol),
#		param_H[:KaneMele_SOC],
#		,
#		tol = hopp_cutoff))

        # -------------- Rashba spin-orbit coupling -------------- #      
#
#  push!(hoppings, TBmodel.Hopping_Term( 
#		(ri,rj) -> isapprox(LA.norm(ri.-rj),dist_nnn,atol=dist_tol),
#		param_H[:Rashba_SOC],
#		,
#		tol = hopp_cutoff))
#


        # -------------- compute max_dist and sum up ------------- #      



  insheet(ri,rj) = (dim==2 || same([ri[3],rj[3]],AtomsUC[1:2,3])) && LA.norm(ri[1:2].-rj[1:2]) < maxdist
  

  nr = count( Algebra.OuterDist(AllAtoms[:,1:2],AtomsUC[1:1,1:2]) .< maxdist)


  return TB_Input, (Sum(Hoppings,insheet), nr, hopp_cutoff)



end

#===========================================================================#
#
# Collect Hamilt terms into one
#
#---------------------------------------------------------------------------#


function OperatorData(Lattice,param_sheet,LattLabels=[1],param_coupling=();dist_tol=1e-5,hopp_cutoff=1e-6)

  sheet(iL,L) = Sheet(param_sheet(iL),L,dist_tol,hopp_cutoff)

  length(LattLabels)==1 && return sheet(1,Lattice)
  


  TB_Input, (h_c, nr_c_, aux) = Coupling_Sheets(param_coupling,Lattice,dist_tol,hopp_cutoff)

  nr_c = nr_c_*(length(param_coupling) > 0)


  sheets = [sheet(L,Lattice.Select_Sublattices(L))[2] for L in LattLabels]

  H =  Utils.Sum_functions(h_c,[h for (h,nr,hc) in sheets]...)

  Nr = length(sheets)*(maximum([nr for (h,nr,hc) in sheets]) + nr_c)


  return (TB_Input,(H,Nr,hopp_cutoff))


end


#===========================================================================#
#
# Valley operator
#
#---------------------------------------------------------------------------#

#data from file or produce it
function Valley_Operator(savefile::String="";argH::String="")

  O = TBmodel.Bloch_Hamilt(savefile,argH=argH)

  return (P;k,kwargs...)->reshape(real(LA.diag(P'*O(k)*P)),:,1)

end



function Valley_Operator(Lattice,LattLabels=[1];dist_tol=1e-5,hopp_cutoff=1e-6,parallel,argH,savefile)


  data = OperatorData(Lattice,L->(Intralayer_Hopping=0.0,AntiHaldane=1.0,),LattLabels,dist_tol=dist_tol,hopp_cutoff=hopp_cutoff)


  O = TBmodel.Bloch_Hamilt(data...,parallel=parallel,argH=argH,savefile=savefile)


  return (P;k,kwargs...)-> reshape(LA.diag(P'*O(k)*P),:,1)

end


#===========================================================================#
#
# Layer operator
#
#---------------------------------------------------------------------------#


function Layer_Operator(Rs::Matrix{Float64},N_Ham::Int64)#,d0=1;savefile="")

  return TBmodel.Position_Operator(3,Rs,N_Ham)
#  Z = Rs[:,3:3]
#
#  f_ldos = TBmodel.LDOS(size(Rs,1),N_Ham)
#
#  return (P;kwargs...) -> f_ldos(P)*Z

#  function get_ldos(P,ldos)
#  
#    !isnothing(ldos) && return ldos
#  
#    return f_ldos(P)
#  
#  end
#
#  return (P;LDOS=nothing,kwargs...) -> get_ldos(P,LDOS)*Z

end




#  f(i) = TBmodel.flat_indsvals(i,i,Matrix(LA.I,d0,d0).*Rs[i,3],d0)
#
#  T0 = TBmodel.indsvals_to_Tm(size(Rs,1)*d0,vcat(f.(axes(Rs,1))...))
#
##  TBmodel.Bloch_Hamilt(data...,kwargs...)
#
#  data = ([ zeros(Int64,size(Rs,2)) ],[ zeros(Float64,size(Rs,2)) ],[T0],)
#
#  TBmodel.Save_Hamilt_Data(savefile,data...)
#
#  return TBmodel.Assemble_Bloch_Hamiltonian(data)





































#===========================================================================#
#
# Finds the Dirac point(s) - doesn't really work
#
#---------------------------------------------------------------------------#

function Dirac_Point(H,K,nr_bands=nothing,tol=1e-8,window=0.05)

  !isnothing(nr_bands)  && (nr_bands *= 3)

  EnK = BandStructure.get_eigen(H,tol=tol,nr_bands=nr_bands)(K).Es

  EnK = EnK[abs.(EnK) .< window]

  D = Algebra.OuterDist(EnK,EnK)

  
	# Catesian Index [CI] containing the two closest eigenvalues,


  indices = findall(LA.triu(trues(size(D)...),1))
#  indices = findall(LA.triu(trues(length(EnK),length(EnK)),1))


#  indices[argmin(D[indices])]

  D =  D[indices] 

  CIs = indices[findall(abs.(D .- minimum(D)) .< tol)]

   
 
#  println(findall(LA.triu(abs.(D.-minimum(D[LA.triu(trues(size(D)...),1)])) .< tol,1)))


#  CIs = findall(LA.triu(abs.(D.-minimum(D[findall(LA.triu(D.>0))])) .< tol,1))

  return [EnK[i] for i in Tuple(CIs[1])] |> E -> real(sum(E)/length(E))

#  if length(CIs) in [1,2]
#
#    return [EnK[i] for p in map(Tuple,CIs) for i in p] |> E -> sum(E)/length(E)
#
#  end

#  for w in range(0,maximum(abs.(EnK)),length=100)
#
#  EnK = EnK[abs.(EnK) .< window]
#
#  d1,d2 = 
#	# Catesian Indices CIs
#
#  for maxendiff in collect(100:-5:1)*minimum()
#
#    CIs = findall(LA.triu(D.<maxendiff,1))
#
#
#    !mute && println(length(CIs), "diffs ",[-(EnK[collect(Tuple(ci))]...) for ci in CIs])
#
#
#  end

  error("The Dirac points couldn't be found/too many found.")




end


###############################################################################

end
