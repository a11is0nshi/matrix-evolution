# Given R, this will return <e
from networkx import *
import matplotlib.pyplot as plt 

txt_file = "AML-10_k_1.txt"
f = open(txt_file).read()
# this will be the set that is read in from the text file
s = eval(f)

R = {(15, 21), (38, 23), (7, 17), (18, 17), (15, 30), (7, 26), (18, 26), (15, 39), (7, 35), (18, 35), (29, 32), (8, 9), (0, 5), (8, 18), (0, 14), (11, 14), (21, 46), (0, 23), (11, 23), (33, 20), (44, 20), (44, 29), (4, 2), (44, 38), (2, 32), (25, 34), (33, 47), (22, 10), (44, 47), (25, 43), (3, 6), (34, 21), (45, 21), (37, 17), (34, 30), (45, 30), (3, 24), (37, 26), (3, 33), (37, 35), (22, 46), (3, 42), (38, 9), (7, 3), (15, 16), (7, 21), (18, 21), (7, 30), (18, 30), (48, 40), (41, 1), (21, 32), (0, 9), (6, 43), (41, 10), (44, 15), (33, 24), (44, 24), (25, 20), (44, 33), (25, 29), (33, 42), (44, 42), (25, 38), (3, 1), (14, 1), (34, 16), (45, 16), (25, 47), (3, 10), (14, 10), (22, 23), (3, 19), (22, 32), (3, 37), (15, 2), (36, 43), (28, 39), (7, 16), (18, 16), (48, 17), (5, 46), (48, 26), (29, 22), (48, 35), (29, 40), (33, 1), (6, 47), (44, 1), (33, 10), (44, 10), (41, 23), (25, 6), (33, 19), (44, 19), (41, 32), (25, 15), (44, 28), (25, 24), (33, 37), (34, 2), (44, 37), (25, 33), (45, 2), (2, 40), (3, 5), (25, 42), (3, 14), (13, 46), (24, 46), (3, 23), (14, 23), (14, 32), (15, 6), (7, 2), (18, 2), (28, 34), (36, 47), (9, 30), (28, 43), (29, 17), (6, 24), (29, 26), (29, 35), (6, 42), (21, 40), (44, 5), (25, 1), (33, 14), (44, 14), (25, 10), (33, 23), (44, 23), (2, 17), (25, 19), (2, 26), (25, 28), (2, 35), (13, 32), (24, 32), (3, 9), (16, 46), (36, 24), (28, 20), (28, 29), (36, 42), (28, 47), (6, 10), (9, 43), (6, 19), (29, 21), (21, 17), (29, 30), (21, 26), (6, 37), (21, 35), (6, 46), (39, 43), (44, 9), (25, 5), (44, 18), (25, 14), (43, 40), (36, 1), (16, 32), (36, 10), (28, 6), (36, 19), (5, 22), (28, 24), (36, 37), (9, 20), (28, 33), (3, 47), (5, 40), (28, 42), (6, 14), (29, 16), (9, 47), (6, 23), (6, 32), (39, 47), (25, 9), (31, 43), (12, 39), (43, 17), (12, 48), (43, 26), (13, 22), (43, 35), (13, 40), (24, 40), (28, 1), (36, 14), (28, 10), (36, 23), (5, 17), (28, 19), (36, 32), (5, 26), (9, 24), (5, 35), (28, 37), (29, 2), (8, 46), (19, 46), (30, 46), (39, 24), (39, 42), (12, 34), (20, 47), (4, 30), (31, 47), (12, 43), (4, 39), (13, 17), (24, 17), (4, 48), (13, 26), (24, 26), (13, 35), (24, 35), (16, 40), (28, 5), (9, 1), (28, 14), (9, 10), (38, 46), (5, 21), (28, 23), (9, 19), (48, 1), (5, 30), (8, 32), (19, 32), (39, 10), (8, 41), (0, 37), (11, 37), (39, 19), (0, 46), (11, 46), (31, 24), (39, 37), (12, 20), (39, 46), (12, 29), (31, 42), (12, 38), (4, 34), (12, 47), (4, 43), (13, 21), (16, 17), (13, 30), (16, 26), (16, 35), (38, 32), (28, 9), (5, 16), (9, 14), (27, 22), (8, 27), (27, 40), (8, 36), (20, 1), (31, 1), (0, 32), (11, 32), (39, 14), (20, 10), (31, 10), (0, 41), (12, 6), (39, 23), (31, 19), (42, 47), (39, 32), (12, 15), (12, 24), (31, 37), (4, 20), (12, 33), (4, 29), (12, 42), (4, 38), (13, 16), (4, 47), (38, 27), (5, 2), (15, 34), (38, 36), (15, 43), (7, 39), (18, 39), (27, 17), (7, 48), (8, 13), (27, 26), (8, 22), (27, 35), (0, 18), (8, 31), (0, 27), (8, 40), (19, 40), (30, 40), (0, 36), (12, 1), (31, 14), (12, 10), (41, 46), (20, 23), (4, 6), (31, 23), (12, 19), (20, 32), (4, 15), (31, 32), (12, 28), (4, 24), (12, 37), (13, 2), (4, 33), (34, 43), (45, 43), (3, 46), (14, 46), (38, 13), (15, 20), (38, 22), (15, 29), (38, 31), (38, 40), (7, 34), (18, 34), (15, 47), (42, 1), (7, 43), (18, 43), (42, 10), (27, 21), (8, 17), (19, 17), (27, 30), (0, 13), (30, 17), (8, 26), (19, 26), (30, 26), (0, 22), (11, 22), (8, 35), (19, 35), (30, 35), (0, 31), (0, 40), (11, 40), (12, 5), (39, 22), (4, 1), (12, 14), (4, 10), (33, 46), (12, 23), (44, 46), (4, 19), (34, 20), (43, 1), (45, 20), (4, 28), (3, 32), (34, 47), (45, 47), (38, 17), (15, 24), (38, 26), (7, 20), (18, 20), (15, 33), (38, 35), (7, 29), (18, 29), (15, 42), (7, 38), (8, 3), (27, 16), (7, 47), (18, 47), (8, 21), (42, 23), (0, 17), (11, 17), (8, 30), (0, 26), (11, 26), (0, 35), (11, 35), (33, 32), (44, 32), (12, 9), (4, 5), (44, 41), (25, 37), (4, 14), (25, 46), (34, 24), (45, 24), (3, 27), (22, 40), (34, 42), (45, 42), (3, 36), (15, 1), (37, 47), (15, 10), (7, 6), (18, 6), (15, 19), (38, 21), (7, 15), (38, 30), (7, 24), (18, 24), (15, 37), (27, 2), (7, 33), (18, 33), (7, 42), (18, 42), (27, 20), (0, 3), (8, 16), (0, 21), (11, 21), (11, 30), (44, 27), (41, 40), (25, 23), (44, 36), (34, 1), (45, 1), (25, 32), (4, 9), (34, 10), (45, 10), (25, 41), (22, 17), (34, 19), (45, 19), (3, 13), (22, 26), (3, 22), (22, 35), (34, 37), (45, 37), (3, 31), (3, 40), (14, 40), (15, 5), (7, 1), (18, 1), (15, 14), (38, 16), (7, 10), (18, 10), (15, 23), (36, 46), (7, 19), (18, 19), (15, 32), (7, 28), (7, 37), (8, 2), (18, 37), (48, 47), (29, 43), (0, 16), (11, 16), (41, 17), (44, 13), (41, 26), (33, 22), (44, 22), (41, 35), (25, 18), (44, 31), (25, 27), (33, 40), (44, 40), (25, 36), (37, 1), (34, 14), (45, 14), (37, 10), (34, 23), (45, 23), (3, 17), (14, 17), (45, 32), (3, 26), (14, 26), (3, 35), (14, 35), (38, 2), (15, 9), (7, 5), (18, 5), (7, 14), (18, 14), (28, 46), (7, 23), (9, 42), (18, 23), (27, 1), (29, 20), (0, 2), (11, 2), (29, 47), (33, 17), (44, 17), (25, 13), (33, 26), (44, 26), (25, 22), (33, 35), (44, 35), (25, 31), (25, 40), (2, 47), (3, 21), (37, 23), (3, 30), (15, 13), (28, 32), (7, 9), (18, 9), (48, 10), (9, 37), (9, 46), (6, 22), (29, 24), (48, 46), (6, 40), (29, 42), (44, 3), (21, 47), (33, 21), (44, 21), (25, 17), (33, 30), (44, 30), (25, 26), (25, 35), (3, 16), (36, 22), (28, 27), (36, 40), (9, 23), (28, 36), (29, 1), (9, 32), (5, 43), (29, 10), (48, 23), (6, 17), (29, 19), (48, 32), (6, 26), (6, 35), (29, 37), (2, 1), (25, 3), (33, 16), (44, 16), (2, 10), (20, 46), (31, 46), (25, 21), (3, 2), (43, 47), (13, 43), (36, 17), (28, 13), (36, 26), (5, 20), (28, 22), (36, 35), (28, 31), (28, 40), (21, 1), (5, 47), (29, 14), (21, 10), (6, 21), (29, 23), (6, 30), (33, 2), (44, 2), (25, 16), (12, 46), (2, 23), (4, 42), (22, 1), (13, 20), (13, 47), (24, 47), (36, 21), (28, 17), (36, 30), (5, 24), (28, 26), (9, 22), (28, 35), (5, 42), (9, 40), (6, 16), (42, 46), (21, 23), (39, 40), (44, 6), (25, 2), (12, 32), (43, 10), (12, 41), (4, 37), (4, 46), (13, 24), (43, 46), (13, 42), (5, 1), (36, 16), (16, 47), (5, 10), (5, 19), (28, 21), (9, 17), (28, 30), (9, 26), (5, 37), (6, 2), (9, 35), (42, 32), (27, 43), (6, 20), (8, 39), (39, 17), (8, 48), (39, 26), (31, 22), (39, 35), (12, 18), (12, 27), (20, 40), (4, 23), (31, 40), (12, 36), (13, 1), (24, 1), (4, 32), (13, 10), (24, 10), (4, 41), (43, 23), (13, 19), (43, 32), (13, 37), (36, 2), (36, 20), (38, 39), (5, 14), (28, 16), (15, 46), (5, 23), (9, 21), (5, 32), (8, 34), (27, 47), (0, 30), (8, 43), (0, 39), (39, 21), (20, 17), (31, 17), (0, 48), (12, 13), (39, 30), (20, 26), (31, 26), (12, 22), (20, 35), (4, 18), (31, 35), (12, 31), (4, 27), (12, 40), (4, 36), (16, 1), (13, 14), (16, 10), (34, 46), (13, 23), (24, 23), (45, 46), (28, 2), (38, 34), (38, 43), (9, 16), (7, 46), (18, 46), (27, 24), (6, 1), (8, 20), (8, 29), (27, 42), (8, 38), (42, 40), (0, 34), (39, 16), (8, 47), (19, 47), (30, 47), (0, 43), (11, 43), (31, 21), (12, 17), (31, 30), (4, 13), (12, 26), (4, 22), (12, 35), (4, 31), (34, 32), (4, 40), (37, 46), (16, 23), (38, 20), (15, 27), (38, 29), (9, 2), (15, 36), (7, 32), (18, 32), (38, 47), (27, 10), (7, 41), (8, 6), (27, 19), (8, 15), (42, 17), (8, 24), (42, 26), (27, 37), (0, 20), (11, 20), (8, 33), (39, 2), (27, 46), (0, 29), (42, 35), (8, 42), (0, 38), (12, 3), (39, 20), (31, 16), (0, 47), (11, 47), (12, 21), (4, 17), (12, 30), (4, 26), (4, 35), (37, 32), (3, 39), (38, 6), (15, 22), (38, 24), (7, 18), (15, 31), (38, 33), (7, 27), (18, 27), (15, 40), (38, 42), (7, 36), (8, 1), (18, 36), (19, 1), (27, 14), (30, 1), (8, 10), (19, 10), (27, 23), (0, 6), (30, 10), (8, 19), (27, 32), (0, 15), (8, 28), (0, 24), (11, 24), (8, 37), (31, 2), (0, 33), (0, 42), (11, 42), (31, 20), (4, 3), (44, 39), (12, 16), (44, 48), (4, 21), (34, 22), (45, 22), (34, 40), (45, 40), (3, 34), (38, 1), (22, 47), (3, 43), (38, 10), (15, 17), (38, 19), (7, 13), (18, 13), (15, 26), (7, 22), (18, 22), (15, 35), (38, 37), (7, 31), (18, 31), (7, 40), (8, 5), (18, 40), (0, 1), (11, 1), (8, 14), (0, 10), (11, 10), (8, 23), (19, 23), (29, 46), (0, 19), (11, 19), (30, 23), (30, 32), (39, 1), (0, 28), (12, 2), (44, 34), (41, 47), (25, 30), (33, 43), (44, 43), (25, 39), (4, 16), (34, 17), (45, 17), (2, 46), (25, 48), (34, 26), (45, 26), (3, 20), (34, 35), (45, 35), (3, 29), (37, 40), (38, 5), (14, 47), (38, 14)}
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

 # ex. (0,1), (1,2), (0,2)  
G = transitive_reduction(G)

plt.figure(figsize=(12, 12))
pos = spring_layout(G)  # positions for all nodes

draw(G, pos, with_labels=True, node_size=500, node_color="skyblue", font_size=10, font_weight="bold", arrows=True)
plt.title("Directed Graph Visualization")
plt.show()
# f = open("AML-10_k_1.txt", "a")
# f.write(str(G.edges()))
# print(f"edges post trans red: {G.edges}")
