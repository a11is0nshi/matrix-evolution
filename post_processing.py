# Given R, this will return <e
from networkx import *
import graphviz

# txt_file = "nodups_AML-10_k_0.txt"
# f = open(txt_file).read()
# this will be the set that is read in from the text file
# s = eval(f)

dataSet = "smalltest"
R = {(1, 0), (1, 2), (0, 3), (1, 3)}

beta = 20
color = 'aquamarine2'
# filename = "edgelist.txt"
# input: set R, output list of edges
def prune(R):
    Rlist = list(R)
    G = DiGraph(Rlist)
    for scc in strongly_connected_components(G):
        l1 = [num for num in scc]
        l2 = [str(scc)[1:len(str(scc))-1] for i in range(len(scc))]
        mapping = dict(zip(l1, l2))
        G = relabel_nodes(G, mapping)
    G.remove_edges_from(selfloop_edges(G))
    print(G.nodes)
    print(G.edges)
    G = transitive_reduction(G)   # ex. (0,1), (1,2), (0,2) 
    return G
  
def genGraph():
    G = prune(R)
    dot = graphviz.Digraph(f"{dataSet}_beta{beta}",
                           
                            node_attr={'color': 'black', 'style': 'filled', 'fillcolor': color, 'shape': 'circle'},
                            edge_attr={'arrowsize': '0.3'})  
    # Graph color namse: https://graphviz.org/doc/info/colors.html
    dot.attr(dpi='1000')
    for node in G.nodes():
        dot.node(node, label=None)
    for edge in G.edges():
        dot.edge(str(edge[0]), str(edge[1]))
    # dot.view()
    # Save and render the graph
    dot.render(f"{dataSet}_beta{beta}", format='png', cleanup=True)


# def edgesToFile():
#   G = prune(R)
#   f = open("postProEdges_" + dataSet + "_k" + {k}, "a")
#   f.write(str(G.edges()))
#   print(f"edges post trans red: {G.edges}")


genGraph()
