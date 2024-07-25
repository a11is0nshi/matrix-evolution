# Given R, this will return <e
from networkx import *
import graphviz

# txt_file = "nodups_AML-10_k_0.txt"
# f = open(txt_file).read()
# this will be the set that is read in from the text file
# s = eval(f)

data_set = "testNodes"
R = {(0, 1), (1,0), (1,2),(0,2),(1,3), (0,3)}
k = 2
# filename = "edgelist.txt"
# input: set R, output list of edges
def prune(R):
  Rlist = list(R)
  G = DiGraph(Rlist)
  G.remove_edges_from(selfloop_edges(G))   # ex. (1,1)
  remove = []
  for scc in strongly_connected_components(G):
      if len(scc) > 1:
          l = list(scc)
          for i in range(len(l)):
              if i != 0:
                  remove.append(l[i])
              if i == 0:
                  print(f"l[i] = {l[i]}")
                  print(f"scc = {scc}")
                  G = relabel_nodes(G, {l[0]: str(scc)})
                  print(f"G.nodes = {G.nodes()}")
                 
  G.remove_nodes_from(remove)     # ex. (1,0), (0,1)
  G = transitive_reduction(G)   # ex. (0,1), (1,2), (0,2) 
  return G
  
def genGraph():
  G = prune(R)
  dot = graphviz.Digraph(filename=f"{data_set}_k-{k}.gv")
  for edge in G.edges():
      dot.edge(str(edge[0]), str(edge[1]))
  # dot.view()
  # Save and render the graph
  dot.render(f"{data_set}_k{k}", format='png', cleanup=True)

"""
def edgesToFile():
  G = prune(R)
  f = open(filename, "a")
  f.write(str(G.edges()))
  print(f"edges post trans red: {G.edges}")
"""

genGraph()
