# Given R, this will return <e
from networkx import *

R = {(0,1),(0,1),(1,1),(0,2),(1,3),(1,4),(2,5),(2,6)}
Rlist = list(R)
G = DiGraph(Rlist)
G.remove_edges_from(selfloop_edges(G)) # ex. (1,1)
TR = transitive_reduction(G) # ex. (0,1), (1,2), (0,2)
print(list(TR.edges))
