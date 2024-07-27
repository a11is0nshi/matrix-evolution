# Given R, this will return <e
from networkx import *
import graphviz

# txt_file = "nodups_AML-10_k_0.txt"
# f = open(txt_file).read()
# this will be the set that is read in from the text file
# s = eval(f)

dataSet = "AML10-7"
R = {(4, 0), (4, 9), (3, 7), (4, 6), (5, 1), (8, 0), (5, 13), (9, 2), (3, 13), (5, 10), (8, 9), (10, 6), (14, 13), (8, 6), (1, 0), (13, 2), (0, 14), (10, 12), (1, 9), (11, 14), (6, 2), (7, 1), (7, 10), (4, 2), (3, 0), (14, 0), (3, 9), (5, 6), (4, 8), (3, 6), (14, 9), (8, 2), (5, 12), (8, 11), (1, 2), (10, 14), (1, 11), (0, 13), (11, 13), (1, 8), (6, 13), (7, 12), (3, 2), (14, 2), (4, 1), (12, 14), (3, 11), (4, 10), (3, 8), (8, 1), (5, 14), (11, 0), (9, 12), (0, 9), (8, 10), (11, 9), (0, 6), (11, 6), (1, 13), (2, 12), (1, 10), (13, 12), (7, 5), (7, 14), (6, 12), (4, 3), (3, 1), (4, 12), (3, 10), (5, 7), (12, 13), (0, 2), (10, 0), (11, 2), (9, 14), (8, 12), (10, 9), (1, 6), (11, 8), (2, 14), (1, 12), (13, 14), (12, 0), (6, 14), (7, 13), (12, 9), (14, 6), (4, 5), (12, 6), (5, 0), (4, 14), (3, 12), (5, 9), (4, 11), (14, 12), (10, 2), (11, 1), (8, 14), (10, 11), (9, 13), (11, 10), (10, 8), (7, 0), (1, 14), (2, 13), (7, 9), (7, 6), (12, 2), (4, 7), (3, 5), (5, 2), (9, 0), (3, 14), (5, 11), (4, 13), (5, 8), (9, 6), (2, 0), (10, 1), (10, 13), (13, 0), (0, 12), (2, 9), (8, 13), (11, 12), (2, 6), (13, 9), (7, 2), (6, 0), (13, 6), (7, 11), (6, 9), (7, 8)}
beta = 3
color = 'plum3'
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
