# Given R, this will return <e
from networkx import *

R = { }
Rlist = list(R)
G = DiGraph(Rlist)

# ex. (1,1)
G.remove_edges_from(selfloop_edges(G)) 

# ex. (1,0), (0,1)
remove = []
for scc in strongly_connected_components(G):
    if len(scc) ==  2:
      remove.append(list(scc)[0])

for n in remove:
   G.remove_node(n)

 # ex. (0,1), (1,2), (0,2)  
G = transitive_reduction(G)

print(G.edges)
