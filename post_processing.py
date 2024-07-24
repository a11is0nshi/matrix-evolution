# Given R, this will return <e
import networkx as nx

R = {(0,1), (0,2),(1,3),(1,4),(2,5),(2,6)}
Rlist = list(R)
G = nx.Graph(Rlist)
