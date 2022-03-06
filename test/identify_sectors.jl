using Distributed 
include("identify_sectors_.jl")


@show length(Bonds) 


import myLibs: Utils 


pmap(1:20) do i

	for n in ([100,1000,10000])

Utils.IdentifySectors(first, Bonds[1:n]);

	end 
		
	Utils.IdentifySectors(first, Bonds)



end 
