include("rgf.jl")


PyPlot.close(11)


fig,Ax = PyPlot.subplots(3,2,num=11,figsize=(10,11),
												 sharey="col")


nr_atoms_layer_ = [12,24,36,72]

ENERGIES = LinRange(0.001,3,100) 


for (i_width,nr_atoms_layer) in enumerate(nr_atoms_layer_)

for i_latt_type  in 1:2 





lab, get_latt = [("square", square_latt_lead),
										 ("honey", honeycomb_latt_lead)
										 ][i_latt_type]
































fig.suptitle(string("Figure ",fig.number))


fig.tight_layout()
fig.subplots_adjust(top=0.95)
