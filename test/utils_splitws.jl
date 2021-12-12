import myLibs: Utils,Algebra
import Random 
const SW = Utils.SplitWorkstations

using OrderedCollections 


function mem_needed(L::Real)::Real
	
	estimated_memory_needed = Dict(40=>1.,
															 50=>1.7,
															 60=>2,
															 70=>3.4,
															 80=>3.8,
															 90=>5.2,
															 100=>7.2,
															 110=>11,
															 120=>22) 

	K = sort(collect(keys(estimated_memory_needed)))

	V = [estimated_memory_needed[k] for k in K]*1.4

	L > K[end] && return V[end]

	return V[findfirst(>=(L), K)]

end 

println()

possible_machines = [OrderedDict("toad" => 48, "yoshi" => 28, "taranis" => 16, "nut" => 12, "shu" => 11, "re" => 11, "sia" => 8, "neper" => 8, "hapi" => 8), OrderedDict("toad" => 48, "yoshi" => 28, "taranis" => 16, "nut" => 12, "shu" => 11, "re" => 11, "sia" => 8, "neper" => 8, "hapi" => 8), OrderedDict("toad" => 48, "yoshi" => 28, "taranis" => 16, "nut" => 12, "shu" => 11, "re" => 11, "sia" => 8, "neper" => 8, "hapi" => 8), OrderedDict("toad" => 48, "yoshi" => 28, "taranis" => 16, "nut" => 12, "shu" => 11, "re" => 11, "sia" => 8, "neper" => 8, "hapi" => 8), OrderedDict("toad" => 48, "yoshi" => 28, "taranis" => 16, "nut" => 12, "shu" => 11, "re" => 11, "sia" => 8, "neper" => 8, "hapi" => 8), OrderedDict("toad" => 48, "yoshi" => 28, "taranis" => 16, "nut" => 12, "shu" => 11, "re" => 11, "sia" => 8, "neper" => 8, "hapi" => 8), OrderedDict("toad" => 48, "yoshi" => 28, "taranis" => 16, "nut" => 12, "shu" => 11, "re" => 11, "sia" => 8, "neper" => 8, "hapi" => 8), OrderedDict("toad" => 48, "yoshi" => 28, "taranis" => 16, "nut" => 12, "shu" => 11, "re" => 11, "sia" => 8, "neper" => 8, "hapi" => 8), OrderedDict("toad" => 48, "yoshi" => 28, "taranis" => 16, "nut" => 12, "shu" => 11, "re" => 11, "sia" => 8, "neper" => 8, "hapi" => 8), OrderedDict("toad" => 48, "yoshi" => 28, "taranis" => 16, "nut" => 12, "shu" => 11, "re" => 11, "sia" => 8, "neper" => 8, "hapi" => 8), OrderedDict("toad" => 48, "yoshi" => 28, "taranis" => 16, "nut" => 12, "shu" => 11, "re" => 11, "sia" => 8, "neper" => 8, "hapi" => 8), OrderedDict("toad" => 48, "yoshi" => 28, "taranis" => 16, "nut" => 12, "shu" => 11, "re" => 11, "sia" => 8, "neper" => 8, "hapi" => 8), OrderedDict("toad" => 48, "yoshi" => 28, "taranis" => 16, "nut" => 12, "shu" => 11, "re" => 11, "sia" => 8, "neper" => 8, "hapi" => 8), OrderedDict("toad" => 48, "yoshi" => 28, "taranis" => 13, "nut" => 12, "shu" => 6, "re" => 6, "sia" => 6, "neper" => 6, "hapi" => 6), OrderedDict("toad" => 48, "yoshi" => 28, "taranis" => 13, "nut" => 12, "shu" => 6, "re" => 6, "sia" => 6, "neper" => 6, "hapi" => 6), OrderedDict("toad" => 48, "yoshi" => 28, "taranis" => 13, "nut" => 12, "shu" => 6, "re" => 6, "sia" => 6, "neper" => 6, "hapi" => 6), OrderedDict("toad" => 48, "yoshi" => 28, "taranis" => 13, "nut" => 12, "shu" => 6, "re" => 6, "sia" => 6, "neper" => 6, "hapi" => 6), OrderedDict("toad" => 48, "yoshi" => 28, "nut" => 12, "taranis" => 11, "shu" => 5, "re" => 5, "sia" => 5, "neper" => 5, "hapi" => 5), OrderedDict("toad" => 48, "yoshi" => 28, "nut" => 12, "taranis" => 11, "shu" => 5, "re" => 5, "sia" => 5, "neper" => 5, "hapi" => 5), OrderedDict("toad" => 48, "yoshi" => 28, "nut" => 12, "taranis" => 11, "shu" => 5, "re" => 5, "sia" => 5, "neper" => 5, "hapi" => 5), OrderedDict("toad" => 48, "yoshi" => 28, "nut" => 12, "taranis" => 11, "shu" => 5, "re" => 5, "sia" => 5, "neper" => 5, "hapi" => 5), OrderedDict("toad" => 48, "yoshi" => 28, "nut" => 12, "taranis" => 6), OrderedDict("toad" => 48, "yoshi" => 28, "nut" => 12, "taranis" => 6), OrderedDict("toad" => 48, "yoshi" => 28, "nut" => 12, "taranis" => 6), OrderedDict("toad" => 48, "yoshi" => 28, "nut" => 12, "taranis" => 6), OrderedDict("toad" => 48, "yoshi" => 28, "nut" => 11, "taranis" => 5), OrderedDict("toad" => 48, "yoshi" => 28, "nut" => 11, "taranis" => 5), OrderedDict("toad" => 48, "yoshi" => 28, "nut" => 11, "taranis" => 5), OrderedDict("toad" => 48, "yoshi" => 28, "nut" => 11, "taranis" => 5), OrderedDict("toad" => 48, "yoshi" => 28, "nut" => 8, "taranis" => 4), OrderedDict("toad" => 48, "yoshi" => 28, "nut" => 8, "taranis" => 4), OrderedDict("toad" => 48, "yoshi" => 28, "nut" => 8, "taranis" => 4), OrderedDict("toad" => 48, "yoshi" => 28, "nut" => 8, "taranis" => 4), OrderedDict("toad" => 48, "yoshi" => 28, "nut" => 6), OrderedDict("toad" => 48, "yoshi" => 28, "nut" => 6), OrderedDict("toad" => 48, "yoshi" => 28, "nut" => 6), OrderedDict("toad" => 48, "yoshi" => 28, "nut" => 6), OrderedDict("toad" => 48, "yoshi" => 22, "nut" => 4), OrderedDict("toad" => 48, "yoshi" => 22, "nut" => 4), OrderedDict("toad" => 48, "yoshi" => 22, "nut" => 4), OrderedDict("toad" => 48, "yoshi" => 22, "nut" => 4), OrderedDict("toad" => 37, "yoshi" => 11), OrderedDict("toad" => 37, "yoshi" => 11), OrderedDict("toad" => 37, "yoshi" => 11), OrderedDict("toad" => 37, "yoshi" => 11), OrderedDict("toad" => 37, "yoshi" => 11), OrderedDict("toad" => 37, "yoshi" => 11), OrderedDict("toad" => 37, "yoshi" => 11), OrderedDict("toad" => 37, "yoshi" => 11)]
M = [1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 2.38, 2.38, 2.38, 2.38, 2.8, 2.8, 2.8, 2.8, 4.76, 4.76, 4.76, 4.76, 5.319999999999999, 5.319999999999999, 5.319999999999999, 5.319999999999999, 7.279999999999999, 7.279999999999999, 7.279999999999999, 7.279999999999999, 10.08, 10.08, 10.08, 10.08, 15.399999999999999, 15.399999999999999, 15.399999999999999, 15.399999999999999, 30.799999999999997, 30.799999999999997, 30.799999999999997, 30.799999999999997, 30.799999999999997, 30.799999999999997, 30.799999999999997, 30.799999999999997]



M = [1.0, 4.5, 8.2] 

possible_machines = [OrderedDict("yoshi" => 27, "taranis" => 16, "shu" => 12, "re" => 11, "sia" => 8, "kis" => 8, "thot" => 8, "qetesh" => 8, "uneg" => 8, "dedun" => 8, "nephtys" => 8, "menhit" => 7), OrderedDict("kis" => 8, "yoshi" => 6), OrderedDict("kis" => 7)]


s = SW.split_in_tiers(M, possible_machines) 

for si in s 
	println.(si)
	println()
end 


Q = [4, 28, 48, 8]
@show Q 
@assert Utils.PropDistributeBallsToBoxes(div(sum(Q), minimum(Q)), Q)==[1, 7, 12, 2]


@show Utils.PropDistributeBallsToBoxes_cumulRanges(div(sum(Q), minimum(Q)), Q)



println()

error()

x= Utils.flatmap(1:15) do k

	map(1:10) do  j

	Random.seed!(100k+j)
	i = sort(Utils.Random_Items(axes(M,1),k))

#@show k i 
s = SW.split_in_tiers(M[i], possible_machines[i])

@assert !isempty(s) i

return length(s)

#println("************** $k $j ",length(s))

for si in s 
#	println.(si)
#	println()
end 

#println()
end 
end 

@show Algebra.Mean(x.<3)





