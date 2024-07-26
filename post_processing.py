# Given R, this will return <e
from networkx import *
import graphviz

# txt_file = "nodups_AML-10_k_0.txt"
# f = open(txt_file).read()
# this will be the set that is read in from the text file
# s = eval(f)

dataSet = "AML"
R = {(15, 21), (16, 20), (7, 17), (18, 17), (21, 16), (5, 1), (8, 0), (17, 3), (19, 0), (5, 10), (8, 9), (17, 12), (5, 19), (8, 18), (17, 21), (9, 17), (0, 14), (11, 14), (15, 5), (7, 1), (18, 1), (15, 14), (7, 10), (18, 10), (5, 3), (8, 2), (9, 1), (5, 12), (8, 11), (19, 11), (9, 10), (5, 21), (9, 19), (15, 7), (7, 3), (1, 17), (10, 20), (6, 13), (15, 16), (7, 12), (18, 3), (18, 12), (7, 21), (18, 21), (8, 4), (9, 3), (3, 17), (5, 14), (11, 0), (12, 20), (9, 12), (17, 16), (9, 21), (15, 0), (1, 10), (15, 18), (16, 17), (20, 17), (3, 10), (5, 7), (9, 5), (14, 19), (5, 16), (9, 14), (1, 3), (15, 2), (1, 12), (15, 11), (16, 10), (18, 7), (1, 21), (2, 20), (7, 16), (18, 16), (20, 10), (5, 0), (3, 12), (9, 7), (3, 21), (5, 18), (9, 16), (8, 20), (10, 17), (16, 3), (15, 13), (2, 13), (16, 12), (2, 22), (16, 21), (20, 3), (20, 12), (6, 22), (5, 2), (20, 21), (21, 20), (9, 0), (12, 17), (5, 11), (8, 13), (1, 7), (15, 6), (2, 6), (8, 22), (7, 2), (1, 16), (18, 2), (12, 10), (9, 2), (3, 16), (5, 13), (8, 6), (10, 3), (8, 15), (10, 12), (0, 11), (10, 21), (1, 18), (2, 17), (12, 3), (14, 0), (15, 20), (20, 16), (22, 13), (5, 6), (12, 21), (1, 2), (2, 1), (8, 17), (17, 20), (2, 10), (1, 20), (22, 6), (15, 22), (7, 18), (21, 17), (14, 11), (8, 1), (8, 10), (2, 3), (5, 20), (8, 19), (1, 13), (9, 18), (2, 12), (10, 16), (2, 21), (21, 10), (7, 20), (18, 20), (12, 16), (8, 3), (8, 12), (1, 6), (9, 11), (5, 22), (8, 21), (9, 20), (21, 3), (15, 17), (7, 13), (18, 13), (21, 12), (7, 22), (18, 22), (8, 5), (5, 15), (8, 14), (19, 14), (9, 13), (2, 7), (9, 22), (0, 19), (2, 16), (11, 19), (15, 1), (15, 10), (7, 6), (13, 22), (15, 19), (18, 6), (17, 10), (8, 7), (9, 6), (3, 20), (5, 17), (8, 16), (9, 15), (13, 6), (2, 18), (15, 3), (15, 12), (1, 22)}
beta = 2.1
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
  dot = graphviz.Digraph(f"{dataSet}_beta{beta}")
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
