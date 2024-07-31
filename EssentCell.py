"""
EssentCell implementation. Inspired by the ILP formulation in Malikic et. al, 867-868. 

Asks the user for a CSV file containing binary matrix data, transforms this data into 
a binary mutation matrix D, and extracts the essential partial order relation. 

D is a binary n x m matrix, where there are n sequenced single cells and m mutations. 

Additional user-provided parameters: kappa, the option to customize the results (essential 
partial order graph and output file or graph only), the color of the graph nodes. 

Outputs the essential partial order graph in .png format and, if specified, a .txt file 
containing the kappa value, the dimensions of the input binary matrix, the number of calls 
to EssILP, the runtime in seconds, the number of nodes in the graph, the poset width of 
the graph, and the relation itself. 

Authors: Allison Shi and Robyn Burger. 
"""
from gurobipy import *
from sys import * 
import numpy as np
import pandas as pd
import time
from networkx import *
import graphviz

# Ask for name of .csv file and desired kappa value
print("Enter name of csv file containing binary matrix data. (Must end in .csv)")
fileName = input("\nCSV file name: ")
outputName = fileName[:-4]
kappa = int(input("\nEnter kappa value: "))

# Determine whether to output graph and output file or graph only
askVerbose = input("\nEnter 'Y' to generate graph and output file. 'N' for graph only: ")
verbose = askVerbose.lower() == "y"

# Initialize variable to track EssILP calls
count = 0

# Transform csv data to a NumPy matrix (may have duplicates)
df = pd.read_csv(fileName)
E = df.to_numpy()

# Remove duplicates from the input matrix
df.drop_duplicates(inplace=True)
D = df.to_numpy(dtype=int)

# Calculate the multiplicities of each row i in D
# Map the rows of D to the number of times they appear in E
M = {}
for i in range(D.shape[0]):
    dups = 0
    for x in range(E.shape[0]):
        if np.all(D[i] == E[x]):
            dups += 1
    M[i] = dups

# Extract the dimensions of D
n = D.shape[0] 
m = D.shape[1]

# Start time for runtime metric
start_time = time.time()

"""
Returns sigma, the minimum number of bit flips required to make D conflict-free. 
"""
def FindOpt():
    try:
        # Silence console output
        env = Env(empty=True)
        env.setParam("OutputFlag", 0)
        env.start()
        model = Model("min_flip_model", env = env)
        model.Params.LogToConsole = 0

        # X will be constrained to be a conflict-free matrix
        X = model.addMVar((n,m), vtype=GRB.BINARY, name="X")

        # Set objective function
        total = sum(sum(M[i]*(1 - D[i, j])*(X[i, j]) + kappa * M[i] * (D[i, j])*(1 - X[i, j]) for j in range(m)) for i in range(n))
        model.setObjective(total, GRB.MINIMIZE)
        
        # Add additional variables to ensure that X is conflict-free
        B01 = model.addMVar((m,m), vtype=GRB.BINARY, name="B01")
        B10 = model.addMVar((m,m), vtype=GRB.BINARY, name="B10")
        B11 = model.addMVar((m,m), vtype=GRB.BINARY, name="B11")

        # Add constraints to enforce that X is conflict-free for each pair of columns (p, q) 
        model.addConstrs(-1 * X[i, p] + X[i, q] <= B01[p,q] for p in range(m) for q in range(p+1, m) for i in range(n)) # (1)
        model.addConstrs(X[i, p] - X[i, q] <= B10[p,q] for p in range(m) for q in range(p+1, m) for i in range(n))  # (2)
        model.addConstrs(X[i, p] + X[i, q] - 1 <= B11[p,q] for p in range(m) for q in range(p+1, m) for i in range(n))  # (3)
        model.addConstrs(B01[p,q] + B10[p,q] + B11[p,q] <= 2 for p in range(m) for q in range(p+1, m))  # (4)
        
        model.optimize()
        sig = model.ObjVal
        return sig

    # This should never happen
    except GurobiError as ex:
        print(f"*********ERROR*********\n{ex}")
        return -1

"""
Returns False if row u is essentially less than some v in Vset, True otherwise. 
"""
def EssILP(u, Vset, sig):
    # Increase EssILP count variable
    global count 
    count += 1
    try:
        # Silence console output
        env = Env(empty=True)
        env.setParam("OutputFlag", 0)
        env.start()
        model = Model("min_flip_model", env = env)
        model.Params.LogToConsole = 0

        # Initialize ILP variables
        X = model.addMVar((n, m), vtype=GRB.BINARY, name="X")
        B01 = model.addMVar((m, m), vtype=GRB.BINARY, name="B01")
        B10 = model.addMVar((m, m), vtype=GRB.BINARY, name="B10")
        B11 = model.addMVar((m, m), vtype=GRB.BINARY, name="B11")

        # Set objective function
        total = sum(sum(M[i]*(1 - D[i, j])*(X[i, j]) + kappa*M[i]*(D[i, j])*(1 - X[i, j]) for j in range(m)) for i in range(n))
        model.setObjective(total, GRB.MINIMIZE)

        # Add constraints to enforce that X is conflict-free for each pair of columns (p, q)
        # Numbers to the right of each constraint correspond to those in the paper 
        model.addConstrs( X[i, q]- X[i, p] <= B01[p,q] for p in range(m) for q in range(p+1, m) for i in range(n))  # (1)
        model.addConstrs(X[i, p] - X[i, q] <= B10[p,q] for p in range(m) for q in range(p+1, m) for i in range(n))  # (2)
        model.addConstrs(X[i, p] + X[i, q] -1  <= B11[p,q] for p in range(m) for q in range(p+1, m) for i in range(n))  # (3)
        model.addConstrs(B01[p,q] + B10[p,q] + B11[p,q] <= 2 for p in range(m) for q in range(p+1, m))  # (4)
        
        # Add additional vars and constraints to test for the essential partial order relation
        nz = len(Vset)
        z = model.addMVar((m,nz), vtype=GRB.BINARY, name="z")
        model.update()
        for i in range(m):
            for v_index in range(nz):
                v = list(Vset)[v_index]
                model.addConstr(X[u, i] - X[v, i] <= z[i, v_index])         # (6)
                model.addConstr(z[i, v_index] <= (X[u, i] - X[v, i] + 1)/2) # (6)
        
        for v in range(nz):
            model.addConstr(sum(z[i, v] for i in range(m) ) >= 1)   # (7)

        model.addConstr(sum(sum(M[i]*(1 - D[i, j])*(X[i, j]) + kappa * M[i] * (D[i, j])*(1 - X[i, j]) for j in range(m)) for i in range(n)) == sig) # (8)

        model.update()
        model.optimize()

        # The model is infeasible, so u is essentially less than some v in Vset
        if model.Status == 3:
            return False
        # The model is feasible, so u is NOT essentially less than all v in Vset
        else:
            return True
    
    # This should never happen
    except GurobiError as ex:
        print(f"*********ERROR*********\n{ex}")

"""
Splits a set V into two (approximately) equally-sized subsets and returns both subsets. 
"""
def Split(V):
    V = list(V)
    i = int(len(V)/2)
    return set(V[:i]), set(V[i:])

"""
Given a row u and set of rows V, GetRelated(u, V, sig) outputs the set of samples v ∈ V 
such that u is essentially less than v.
"""
def GetRelated(u, V, sig):
    if not EssILP(u, V, sig): 
        if len(V) == 1:
            return V
        else: 
            Vl, Vr = Split(V)
            return GetRelated(u, Vl, sig).union(GetRelated(u, Vr, sig))
    # Returns the empty set if u is NOT essentially less than all v ∈ V
    else: 
        return set()
        
"""
Returns the essential order relation of an input D matrix.
"""
def GetEssential():
    sig = FindOpt()
    # The set of all 0-based row indices
    S = {num for num in range(D.shape[0])}
    # Initialize set of essential order relations
    ess_set = set()

    # Find all essential order relations for each row u
    R = []
    for u in S:
        V = S.difference({u})
        R.append(GetRelated(u, V, sig))
        P = {(u, y) for y in R[u]}
        temp = ess_set.union(P)
        ess_set = temp
    
    # Enforce antisymmetry of the relation by removing all strongly connected components (SCCs) 
    # and replacing each SCC with one node that contains the names of all the nodes in that SCC
    ess_list = list(ess_set)
    G = DiGraph(ess_list)
    for scc in strongly_connected_components(G):
        l1 = [num for num in scc]
        l2 = [str(scc)[1:len(str(scc))-1] for i in range(len(scc))]
        mapping = dict(zip(l1, l2))
        G = relabel_nodes(G, mapping)
    
    # Enforce that there are no self-loop (reflexive) edges in the graph
    G.remove_edges_from(selfloop_edges(G))
    
    # Enforce transitive property of the partial order relation
    G = transitive_reduction(G)
    return G

"""
Generates the essential partial order graph of a specified graphColor. 
"""
def GenGraph(graphColor):
    Graph = GetEssential()
    dot = graphviz.Digraph(f"{outputName}_kappa{kappa}",
                            node_attr={'color': 'black', 'style': 'filled', 'fillcolor': graphColor, 'shape': 'circle'},
                            edge_attr={'arrowsize': '0.3'})  
    # Specify the resolution of the graph 
    dot.attr(dpi='1000')
    
    # Label nodes and edges
    for node in Graph.nodes():
        dot.node(node, label=None)
    for edge in Graph.edges():
        dot.edge(str(edge[0]), str(edge[1]))

    # Save and render the graph
    dot.render(f"{outputName}_kappa{kappa}", format='png', cleanup=True)
    return Graph

"""
Returns poset width of a given graph by calculating the length of the longest antichain. 
"""
def Width(G):
    max = 0
    for x in antichains(G):
        if len(x) > max:
            max = len(x)
    return max

# Name the output results file
fileResults = outputName + "_kappa" + str(kappa) + ".txt"

# Prompt user for the color of the graph nodes
print("\nThe color of the graph nodes must be one of the options listed at \nhttps://graphviz.org/doc/info/colors.html.")
graphColor = input("\nEnter color of graph nodes: ")
G = GenGraph(graphColor)

# Write results to file if user requests file output
if verbose:
    f = open(fileResults, "a")
    f.write(f"kappa: {kappa}\n")
    f.write(f"n (number of samples): {n}\n")
    f.write(f"m (number of mutations): {m}\n")
    f.write(f"EssILP calls: {count}\n")
    f.write(f"Runtime: {time.time() - start_time} seconds\n")
    f.write(f"Number of Nodes: {G.number_of_nodes()}\n")
    f.write(f"Poset Width: {Width(G)}\n")
    f.write(f"Essential Relation: {[edge for edge in G.edges]}\n\n")
    # Inform the user of the name of the results .txt file
    print(f"\nResults were written to: {fileResults}")
    f.close()

# Inform the user of the name of the graph .png file
print(f"\nGraph was generated: {outputName}_kappa{kappa}.png")
