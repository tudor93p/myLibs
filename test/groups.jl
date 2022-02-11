

import myLibs: Groups 


for a in (rand(4),rand(3),
					(rand(4),),(rand(3),)
					)

	Groups.WeylRepr(a...)
	Groups.dWeylRepr(a...)


end 
					
for a in [(rand(), rand(3)), (rand(), rand(1:3)), rand(3), (rand(3),)]


Groups.SU2Repr(a...)
Groups.dSU2Repr(a...)

end 


