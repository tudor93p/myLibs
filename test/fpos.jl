#import myPlots: Transforms 
import myLibs: Operators 

function decode_fstr_Cartesian(s::AbstractString)::Tuple{Int,String}

	for (i,a_)=enumerate('x':'z'), a=(a_,uppercase(a_))

		s=="$a" && return i,"R" 

		for p in [2,3]
				s=="$a^$p" && return i,"R^$p"
		end 
				
		s=="|$a|" && return i,"abs(R)"
		s=="|$a|^2" && return i,"abs2(R)"
		s=="|$a|^3" && return i,"abs(R)^3"

	end  

	return 0,""

end 


function parse_fstr_Cartesian(s::Union{Char,AbstractString}
														 )::Tuple{Int,Function}

	i,f = decode_fstr_Cartesian(string(s))

	return i, (i==0 ? identity : fstr_to_f(f))

end 


function fstr_to_f(s::AbstractString)::Function 

	E = Meta.parse(s)

	return @eval begin 
		function f(R::Real)::Float64 

			$(E)

		end 
	end 

end 



function foo()

	nr_at = 13 
	nr_orb = 4 



	atoms = rand(2,nr_at)

n = "x"

d,f = parse_fstr_Cartesian(n)



kwargs = (nr_at=nr_at, nr_orb=nr_orb, dim=2)

R = Operators.Position(d, atoms; kwargs..., fpos=f) 

#@show 
R(rand(nr_at*nr_orb,2)) 



ipr = Operators.IPR(;kwargs...)

@show ipr(rand(nr_at*nr_orb,7))

end 


foo()
