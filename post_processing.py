# Given R, this will return <e
from networkx import *
import graphviz

# txt_file = "nodups_AML-10_k_0.txt"
# f = open(txt_file).read()
# this will be the set that is read in from the text file
# s = eval(f)

dataSet = "AML-10_rep10"
R = {(12, 4), (4, 0), (4, 9), (5, 1), (8, 0), (3, 13), (5, 10), (8, 9), (10, 6), (9, 8), (0, 5), (11, 5), (14, 13), (0, 14), (2, 11), (11, 14), (13, 8), (6, 2), (7, 1), (6, 11), (7, 10), (4, 2), (3, 6), (5, 3), (8, 2), (9, 1), (5, 12), (8, 11), (9, 10), (0, 7), (2, 4), (11, 7), (13, 1), (1, 8), (13, 10), (6, 4), (7, 3), (6, 13), (7, 12), (3, 8), (14, 8), (8, 4), (9, 3), (11, 0), (5, 14), (9, 12), (0, 9), (11, 9), (13, 3), (1, 10), (13, 12), (7, 5), (7, 14), (3, 1), (14, 1), (12, 13), (3, 10), (5, 7), (14, 10), (9, 5), (0, 2), (11, 2), (9, 14), (1, 3), (13, 5), (1, 12), (13, 14), (12, 6), (14, 3), (5, 0), (3, 12), (5, 9), (4, 11), (14, 12), (9, 7), (11, 4), (10, 8), (1, 5), (13, 7), (6, 1), (7, 0), (1, 14), (2, 13), (7, 9), (12, 8), (3, 5), (5, 2), (14, 5), (9, 0), (3, 14), (5, 11), (4, 13), (10, 1), (13, 0), (8, 13), (1, 7), (13, 9), (2, 6), (7, 2), (7, 11), (12, 1), (12, 10), (3, 7), (5, 4), (4, 6), (14, 7), (9, 2), (5, 13), (8, 6), (10, 3), (1, 0), (13, 2), (10, 12), (1, 9), (0, 11), (2, 8), (13, 11), (7, 4), (6, 8), (12, 3), (3, 0), (14, 0), (3, 9), (5, 6), (4, 8), (14, 9), (10, 5), (1, 2), (0, 4), (2, 1), (13, 4), (10, 14), (1, 11), (0, 13), (2, 10), (11, 13), (6, 10), (12, 5), (3, 2), (14, 2), (4, 1), (12, 14), (3, 11), (14, 11), (4, 10), (8, 1), (8, 10), (10, 7), (1, 4), (0, 6), (2, 3), (11, 6), (1, 13), (2, 12), (6, 3), (6, 12), (12, 7), (3, 4), (14, 4), (4, 3), (4, 12), (8, 3), (10, 0), (8, 12), (10, 9), (1, 6), (0, 8), (2, 5), (9, 11), (11, 8), (2, 14), (6, 5), (12, 0), (6, 14), (7, 13), (12, 9), (14, 6), (4, 5), (4, 14), (8, 5), (10, 2), (9, 4), (0, 1), (11, 1), (8, 14), (10, 11), (9, 13), (0, 10), (2, 7), (11, 10), (6, 7), (7, 6), (12, 2), (12, 11), (4, 7), (5, 8), (8, 7), (10, 4), (9, 6), (0, 3), (2, 0), (11, 3), (10, 13), (0, 12), (2, 9), (11, 12), (13, 6), (6, 0), (6, 9), (7, 8)}
k = 1
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
  dot = graphviz.Digraph(f"{dataSet}_k{k}")
  for edge in G.edges():
      dot.edge(str(edge[0]), str(edge[1]))
  # dot.view()
  # Save and render the graph
  dot.render(f"{dataSet}_k{k}", format='png', cleanup=True)


# def edgesToFile():
#   G = prune(R)
#   f = open("postProEdges_" + dataSet + "_k" + {k}, "a")
#   f.write(str(G.edges()))
#   print(f"edges post trans red: {G.edges}")


genGraph()
