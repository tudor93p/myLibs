module Graph
#############################################################################

import LightGraphs, MetaGraphs, TikzGraphs,TikzPictures



#===========================================================================#
#
#	Bring together functions from LightGraphs and MetaGraph
#
#---------------------------------------------------------------------------#


MetaGraph = MetaGraphs.MetaGraph


filter_vertices = MetaGraphs.filter_vertices

filter_edges = MetaGraphs.filter_edges

function MetaDiPath(n)

	return MetaGraphs.MetaDiGraph(LightGraphs.path_digraph(n))

end


set_prop! = MetaGraphs.set_prop!

set_props! =	MetaGraphs.set_props!

add_vertex! = LightGraphs.add_vertex!
			
add_edge! = LightGraphs.add_edge!

nv = LightGraphs.nv
ne = LightGraphs.ne

edges = LightGraphs.edges 
vertices = LightGraphs.vertices


all_neighbors = LightGraphs.all_neighbors



neighbors =  LightGraphs.neighbors
in_neighbors =  LightGraphs.inneighbors
out_neighbors =  LightGraphs.outneighbors

get_prop = MetaGraphs.get_prop
has_prop = MetaGraphs.has_prop


function node_by_prop!(g,key)

	MetaGraphs.set_indexing_prop!(g,key)

	node(val) = g[val,key]

	node(val...) = node(val)

	return node

end


function get_edges(a::T) where {T}

	T <: LightGraphs.AbstractGraph && return LightGraphs.edges(a)

	typeof(collect(a)[1]) <: LightGraphs.SimpleGraphs.SimpleEdge && return a


end


function SimpleDiGraph(a)

	return LightGraphs.SimpleDiGraphFromIterator(get_edges(a))

end

function SimpleGraph(a)

	return LightGraphs.SimpleGraphFromIterator(get_edges(a))
	
end


function FindPath(g,n1,n2)

	return LightGraphs.a_star(g,n1,n2)

end






#===========================================================================#
#
#	Plot a graph 
#
#---------------------------------------------------------------------------#

function Plot_Graph(g;fname="graph.pdf",nodelabel=string,edgelabel=e->"",colorrule=i->1)

	isempty(fname) && return

	bgcolors = ["green","red","olive","blue","yellow","magenta","brown","cyan","orange"]



	
	color(i) = colorrule(i)==0 ? "white" : bgcolors[colorrule(i)%length(bgcolors)+1]
	
	p = TikzGraphs.plot(
				g,
				nodelabel.(1:nv(g)),
				edge_labels=Dict(e=>edgelabel(e) for e in edges(g)),
				node_style="draw, rounded corners",
				node_styles=Dict(i=>string("fill=",color(i),"!50") for i in 1:nv(g))
							 )
	
	TikzPictures.save(TikzPictures.PDF(fname), p)

end

#===========================================================================#
#
# Lattice to graph
#
#---------------------------------------------------------------------------#


function PlotAtoms_asGraph(Atoms,isbond;kwargs...)

	Plot_Graph(LightGraphs.SimpleGraph(isbond(Atoms,Atoms));kwargs...)
	

end








#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#



































#############################################################################
end

