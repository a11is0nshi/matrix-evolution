# Given R, this will return <e
from networkx import *

R = {(0,1),(1,0),(2,1),(2,0),(5,6),(6,5),(5,7),(7,5),(6,7),(7,6)}
Rlist = list(R)
G = DiGraph(Rlist)

# ex. (1,1)
G.remove_edges_from(selfloop_edges(G)) 

# ex. (1,0), (0,1)
remove = []
for scc in strongly_connected_components(G):
    if len(scc) > 1:
        l = list(scc)
        for i in range(len(l)):
            if i != 0:
                remove.append(l[i])

# remove nodes from the strongly connected components
G.remove_nodes_from(remove)
print(f"nodes: {G.nodes}")

 # ex. (0,1), (1,2), (0,2)  
G = transitive_reduction(G)

print(f"edges: {G.edges}")
