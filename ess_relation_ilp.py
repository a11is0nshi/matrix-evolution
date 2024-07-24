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
# Change name and k to change file
name = "smallest_test.csv"
k = 0 

# Removes duplicate rows 
df = pd.read_csv(name)
df.drop_duplicates(inplace=True)
D = df.to_numpy()
print(str(D))

n = D.shape[0]  # samples/rows
m = D.shape[1]  # mutations/cols

def GetSigma():
    try:
        env = Env(empty=True) # when set to True, silences console outputs
        env.setParam("OutputFlag",0) # when set to 0, silences console outputs
        env.start()
        model = Model("min_flip_model", env = env)
        model.Params.LogToConsole = 0 # when set to 0, silences console outputs

        # X is a conflict free matrix by constraints
        X = model.addMVar((n,m), vtype=GRB.BINARY, name="X")

        # Objective function
        total = sum(sum((1 - D[i, j])*(X[i, j]) + (D[i, j])*(1 - X[i, j]) for j in range(m)) for i in range(n))
        model.setObjective(total, GRB.MINIMIZE)
        model.addConstr(sum(D[i, j] * (1 - X[i, j]) for i in range(n) for j in range(m)) <= k)
        
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
def TestILP(u, Vset, sig):
    try:
        env = Env(empty=True) # when set to True, silences console outputs
        env.setParam("OutputFlag",0) # when set to 0, silences console outputs
        env.start()
        model = Model("min_flip_model", env = env)
        model.Params.LogToConsole = 0 # when set to 0, silences console outputs
    
        X = model.addMVar((n, m), vtype=GRB.BINARY, name="X") # Conflict free matrix
        B01 = model.addMVar((m, m), vtype=GRB.BINARY, name="B01")
        B10 = model.addMVar((m, m), vtype=GRB.BINARY, name="B10")
        B11 = model.addMVar((m, m), vtype=GRB.BINARY, name="B11")

        total = sum(sum((1 - D[i, j])*(X[i, j]) + (D[i, j])*(1 - X[i, j]) for j in range(m)) for i in range(n))
        model.setObjective(total, GRB.MINIMIZE)

        # Enforce no conflicts by checking each pair of columns (p, q)) 
        model.addConstrs( X[i, q]- X[i, p] <= B01[p,q] for p in range(m) for q in range(p+1, m) for i in range(n)) # (1)
        model.addConstrs(X[i, p] - X[i, q] <= B10[p,q] for p in range(m) for q in range(p+1, m) for i in range(n)) # (2)
        model.addConstrs(X[i, p] + X[i, q] -1  <= B11[p,q] for p in range(m) for q in range(p+1, m) for i in range(n)) # (3)
        model.addConstrs(B01[p,q] + B10[p,q] + B11[p,q] <= 2 for p in range(m) for q in range(p+1, m)) # (4)
        
        # Ensure that there are at most k 1->0 flips
        model.addConstr(sum(D[i, j] * (1 - X[i, j]) for i in range(n) for j in range(m)) <= k) # (5)

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
        model.addConstr(sum(X[i, j] * (1 - D[i, j]) for i in range(n) for j in range(m)) == sig) # (9)

        model.update()
        model.optimize()

        if model.Status == 3:
            print(f"u: {u}, V: {Vset}, sig: {sig} False")
            return False # u <e v
        else:
            print(f"u: {u}, V: {Vset}, sig: {sig} True")
            for var in model.getVars():
                if (var.VarName)[0] == "X":
                    print(f"{var.VarName} = {var.x}")
            return True
 
    except GurobiError as ex:
        print('*********ERROR*********')
        print(f"error: {ex}")


def Split(V):
    V = list(V)
    i = int(len(V)/2)
    return set(V[:i]), set(V[i:])

# Given a sample index u and set of test sample indices V, GetRelated(u, V) 
# outputs all samples v ∈ V such that u <e v
def GetRelated(u, V, sig):
    if not TestILP(u, V, sig): 
        if len(V) == 1:
            return V
        else: 
            Vl, Vr = Split(V)
            return GetRelated(u, Vl, sig).union(GetRelated(u, Vr, sig))
    else: 
        return set()
        
# Given a matrix D and positive integer k, GetEssential(D, k) outputs <e
def GetEssential():
    sig = GetSigma()
    S = {num for num in range(D.shape[0])}
    ess_set = set()
    R = []
    for u in S:
        V = S.difference({u})
        R.append(GetRelated(u, V, sig))
        P = {(u, y) for y in R[u]}
        temp = ess_set.union(P)
        ess_set = temp
    return ess_set
  


print(f"R: {GetEssential()}")

