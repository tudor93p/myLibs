module Graph
#############################################################################

import Graphs, MetaGraphs, TikzGraphs,TikzPictures
import ..Utils



#===========================================================================#
#
#	Bring together functions from Graphs and MetaGraph
#
#---------------------------------------------------------------------------#


MetaGraph = MetaGraphs.MetaGraph


filter_vertices = MetaGraphs.filter_vertices

filter_edges = MetaGraphs.filter_edges

function MetaDiPath(n::Int)

	MetaGraphs.MetaDiGraph(Graphs.path_digraph(n))

end


set_prop! = MetaGraphs.set_prop!

set_props! =	MetaGraphs.set_props!

add_vertex! = Graphs.add_vertex!
			
add_edge! = Graphs.add_edge!

nv = Graphs.nv
ne = Graphs.ne

edges = Graphs.edges 
vertices = Graphs.vertices


all_neighbors = Graphs.all_neighbors



neighbors =  Graphs.neighbors
in_neighbors =  Graphs.inneighbors
out_neighbors =  Graphs.outneighbors

get_prop = MetaGraphs.get_prop
has_prop = MetaGraphs.has_prop


function node_by_prop!(g,key)

	MetaGraphs.set_indexing_prop!(g,key)

	node(val) = g[val,key]

	node(val...) = node(val)

	return node

end


function get_edges(a::Graphs.AbstractGraph)

	Graphs.edges(a)

end 


function get_edges(a::T)::T where T<:Graphs.AbstractEdgeIter

	a

end 

function get_edges(a::Utils.List)::Vector{Graphs.Edge}

	@assert Utils.isList(a, Graphs.AbstractEdge) 

	return collect(a)

end




function SimpleDiGraph(a)

	Graphs.SimpleDiGraphFromIterator(get_edges(a))

end

function SimpleGraph(a)

	Graphs.SimpleGraphFromIterator(get_edges(a))
	
end


function FindPath(g,n1,n2)

	Graphs.a_star(g,n1,n2)

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

	Plot_Graph(Graphs.SimpleGraph(isbond(Atoms,Atoms));kwargs...)
	

end








#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#



































#############################################################################
end

