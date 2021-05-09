using myLibs: TBmodel 
using LinearAlgebra 


R = rand(2,3);
R= hcat(R,R.+[0,1],R.+[1,0])#,R.+[0,-1],R.+[-1,0])

hopp(ri,rj) = isapprox(norm(ri-rj),0)*2.0+isapprox(norm(ri-rj),1)*1.0




#println(TBmodel.HoppingMatrix(R; Hopping=hopp))

TBL = ([0,1],[0,[2,0]],R)

inter, (ms,Rms,Tms) = TBmodel.Compute_Hopping_Matrices(TBL,Hopping=hopp)


@show size(inter)

@show ms 
@show Rms 

@show size.(Tms)



H = TBmodel.Bloch_Hamilt(TBL; Hopping=hopp)


@show H(1)-H([2])+H([3,4])-H([5,6,7])










nothing
