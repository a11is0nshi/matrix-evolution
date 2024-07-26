"""
ILP Implementation inspired by (Malikic et. al, 867-868) paper. 
This implementation sets D to be a random binary m x n matrix, where there are 
m sequenced single cells and n mutations
"""
from gurobipy import *
from sys import * 
import numpy as np
import pandas as pd
import csv
import time as t
from networkx import *
import graphviz

patient = "AML10_"
beta = 1

filename = patient + "B" + str(beta) +  "_results.txt"
f = open(filename, "w")
f.write(f"Name \t Calls \t Width \t Nodes \t Sigma \t time \t Edges\n\n")

def GetSigma(D, M):
    try:
        n = D.shape[0]  # samples/rows
        m = D.shape[1]  # mutations/cols
        env = Env(empty=True) # when set to True, silences console outputs
        env.setParam("OutputFlag",0) # when set to 0, silences console outputs
        env.start()
        model = Model("min_flip_model", env = env)
        model.Params.LogToConsole = 0 # when set to 0, silences console outputs

        # X is a conflict free matrix by constraints
        X = model.addMVar((n,m), vtype=GRB.BINARY, name="X")

        # Objective function
        total = sum(sum(M[i]*(1 - D[i, j])*(X[i, j]) + beta * M[i] * (D[i, j])*(1 - X[i, j]) for j in range(m)) for i in range(n))
        model.setObjective(total, GRB.MINIMIZE)
          
        B01 = model.addMVar((m,m), vtype=GRB.BINARY, name="B01")
        B10 = model.addMVar((m,m), vtype=GRB.BINARY, name="B10")
        B11 = model.addMVar((m,m), vtype=GRB.BINARY, name="B11")

        # Constraints (ensures no conflicts by checking each pair of columns (p, q)) 
        model.addConstrs(-1 * X[i, p] + X[i, q] <= B01[p,q] for p in range(m) for q in range(p+1, m) for i in range(n))
        model.addConstrs(X[i, p] - X[i, q] <= B10[p,q] for p in range(m) for q in range(p+1, m) for i in range(n))
        model.addConstrs(X[i, p] + X[i, q] - 1 <= B11[p,q] for p in range(m) for q in range(p+1, m) for i in range(n))
        model.addConstrs(B01[p,q] + B10[p,q] + B11[p,q] <= 2 for p in range(m) for q in range(p+1, m))
          
        model.optimize()
        sig = model.ObjVal
        #print(f"sigma: {sig}")
        return sig

    except GurobiError as ex:
        print('*********ERROR*********')
        print(ex)
        return -1

# Returns True if u<e v and False otherwise
def TestILP(D, M, u, Vset, sig, count):
    #   global count 
    #  if count % 50 == 0:
    #  print(count)
    count = count + 1
    try:
        n = D.shape[0]  # samples/rows
        m = D.shape[1]  # mutations/cols
        env = Env(empty=True) # when set to True, silences console outputs
        env.setParam("OutputFlag",0) # when set to 0, silences console outputs
        env.start()
        model = Model("min_flip_model", env = env)
        model.Params.LogToConsole = 0 # when set to 0, silences console outputs
      
        X = model.addMVar((n, m), vtype=GRB.BINARY, name="X") # Conflict free matrix
        B01 = model.addMVar((m, m), vtype=GRB.BINARY, name="B01")
        B10 = model.addMVar((m, m), vtype=GRB.BINARY, name="B10")
        B11 = model.addMVar((m, m), vtype=GRB.BINARY, name="B11")

        total = sum(sum(M[i]*(1 - D[i, j])*(X[i, j]) + beta*M[i]*(D[i, j])*(1 - X[i, j]) for j in range(m)) for i in range(n))
        model.setObjective(total, GRB.MINIMIZE)

        # Enforce no conflicts by checking each pair of columns (p, q)) 
        model.addConstrs( X[i, q]- X[i, p] <= B01[p,q] for p in range(m) for q in range(p+1, m) for i in range(n)) # (1)
        model.addConstrs(X[i, p] - X[i, q] <= B10[p,q] for p in range(m) for q in range(p+1, m) for i in range(n)) # (2)
        model.addConstrs(X[i, p] + X[i, q] -1  <= B11[p,q] for p in range(m) for q in range(p+1, m) for i in range(n)) # (3)
        model.addConstrs(B01[p,q] + B10[p,q] + B11[p,q] <= 2 for p in range(m) for q in range(p+1, m)) # (4)
          
        # New constraints 
        nz = len(Vset)
        z = model.addMVar((m,nz), vtype=GRB.BINARY, name="z")
        model.update()
        for i in range(m):
            for v_index in range(nz):
                v = list(Vset)[v_index]
                model.addConstr(X[u, i] - X[v, i] <= z[i, v_index])         # (7)
                model.addConstr(z[i, v_index] <= (X[u, i] - X[v, i] + 1)/2) # (7)
          
        for v in range(nz):
            model.addConstr(sum(z[i, v] for i in range(m) ) >= 1) # (8)

        # is it possible to have at most sig 0 -> 1 flips? 
        model.addConstr(sum(sum(M[i]*(1 - D[i, j])*(X[i, j]) + beta * M[i] * (D[i, j])*(1 - X[i, j]) for j in range(m)) for i in range(n)) == sig) # (9)

        model.update()
        model.optimize()

        if model.Status == 3:
            return False, count # u <e v
        else:
            return True, count
  
    except GurobiError as ex:
        print('*********ERROR*********')
        print(f"error: {ex}")

def Split(V):
    V = list(V)
    i = int(len(V)/2)
    return set(V[:i]), set(V[i:])

# Given a sample index u and set of test sample indices V, GetRelated(u, V) 
# outputs all samples v ∈ V such that u <e v
def GetRelated(D, M, u, V, sig, count):
    val, count = TestILP(D, M, u, V, sig, count)
    if not val: 
        if len(V) == 1:
            return V, count
        else: 
            Vl, Vr = Split(V)
            Sl, cl = GetRelated(D, M, u, Vl, sig, 0)
            Sr, cr = GetRelated(D, M, u, Vr, sig, 0)
            new_set = Sl.union(Sr)
            count = count + cl + cr
            return new_set, count
    else: 
        return set(), count
          
# Given a matrix D, GetEssential(D) outputs <e
def GetEssential(D, M):
    sig = GetSigma(D, M)
    S = {num for num in range(D.shape[0])}
    ess_set = set()
    R = []
    start_time = t.time()
    count = 0
    for u in S:
        V = S.difference({u})
        new_set, new_count = GetRelated(D, M, u, V, sig, count)
        count = new_count
        R.append(new_set)
        P = {(u, y) for y in R[u]}
        temp = ess_set.union(P)
        ess_set = temp
    time = t.time() - start_time
    return ess_set, sig, time, count

def prune(R):
    Rlist = list(R)
    G = DiGraph(Rlist)
    for scc in strongly_connected_components(G):
        l1 = [num for num in scc]
        l2 = [str(scc)[1:len(str(scc))-1] for i in range(len(scc))]
        mapping = dict(zip(l1, l2))
        G = relabel_nodes(G, mapping)
    G.remove_edges_from(selfloop_edges(G))
    return transitive_reduction(G) 

def width(G):
    max = 0
    for x in antichains(G):
        if len(x) > max:
            max = len(x)
    return max

# each file has to be named like this: AML10_1.csv
def ess(f, name, beta):
    fileName = name + ".csv"
    df = pd.read_csv(fileName)
    E = df.to_numpy() # has duplicates
    df.drop_duplicates(inplace=True)
    D = df.to_numpy(dtype=int) # no duplicates
    M = {} # maps rows of D to number of times they appear in E

    #Calculate the number of duplicates
    for i in range(D.shape[0]):
        dups = 0
        for x in range(E.shape[0]):
            if np.all(D[i] == E[x]):
                dups += 1
        M[i] = dups
  
    R, sigma, time, count = GetEssential(D, M)
  
    G = prune(R)
    w = str(width(G))
    numNodes = str(G.number_of_nodes())
    edges = str(G.edges)
    f.write(f"{name} \t {count} \t {w} \t {numNodes} \t {sigma} \t {time} \t {edges} \n")
    print(f"completed {name}!")
  
  # print(f"R: {GetEssential()}")
  # print(f"TestILP was called {count} times")

for i in range(1, 2):
    name = patient + str(i)
    ess(f, name, beta)
print("Done!")