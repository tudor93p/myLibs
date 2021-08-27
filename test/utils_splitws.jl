import myLibs.Utils: SplitWorkstations
const SW = SplitWorkstations

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


M = [1.4, 1.4, 2.8, 4.76, 7.279999999999999, 15.399999999999999, 30.799999999999997]
possible_machines = [OrderedDict("toad" => 48, "yoshi" => 28, "taranis" => 16, "nut" => 12, "shu" => 11, "re" => 11, "sia" => 8, "neper" => 8, "hapi" => 8), 
										 OrderedDict("toad" => 48, "yoshi" => 28, "taranis" => 16, "nut" => 12, "shu" => 11, "re" => 11, "sia" => 8, "neper" => 8, "hapi" => 8), 
										 OrderedDict("toad" => 48, "yoshi" => 28, "nut" => 12, "taranis" => 11, "shu" => 5, "re" => 5, "sia" => 5, "neper" => 5, "hapi" => 5), 
										 OrderedDict("toad" => 48, "yoshi" => 28, "nut" => 12, "taranis" => 6), 
										 OrderedDict("toad" => 48, "yoshi" => 28, "nut" => 8, "taranis" => 4), 
										 OrderedDict("toad" => 48, "yoshi" => 22, "nut" => 4), 
										 OrderedDict("toad" => 37, "yoshi" => 11)]


@time s = SW.split_in_tiers(M, possible_machines) 
@time s = SW.split_in_tiers(M, possible_machines)  

for si in s 
	println.(si)
	println()
end 







