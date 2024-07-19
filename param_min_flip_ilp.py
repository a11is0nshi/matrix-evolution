"""
ILP Implementation inspired by (Malikic et. al, 867-868) paper. 
This implementation sets D to be a random binary m x n matrix, where there are 
m sequenced single cells and n mutations

"""
from gurobipy import *
from sys import * 
import numpy as np
import pandas as pd

# Change name var to change file
name = "Patient2_phyolin.csv"

def getMatrix(name):
    df = pd.read_csv(name)
    return df.to_numpy()

# D is input binary matrix
D = getMatrix(name)

def ILPincreased(u, Vset):
    N = len(Vset)   
    if u <= N:
           # num of samples - rows
        M = D.shape[1]  # num of mutations - cols

        # Populate VMatrix - extracts certain rows of D matrix
    
        V = np.zeros((N, M), dtype=int)
        row_index = 0
        for i in Vset:
            V[row_index] = D[i-1]
            row_index = row_index + 1

        try:
            ma = Model("min_flip_model1")
            mb = Model("min_flip_model2")
        
            # X is a conflict free matrix by constraints
            X = ma.addMVar((N, M), vtype=GRB.BINARY, name="X")
            Y = mb.addMVar((N, M), vtype=GRB.BINARY, name="Y")

            # Objective function - row
            totala = sum(sum((1 - V[i, j])*(X[i, j]) for j in range(M)) for i in range(N))
            ma.setObjective(totala, GRB.MINIMIZE)

            # Objective function - row prime
            totalb = sum(sum((1 - V[i, j])*(Y[i, j]) for j in range(M)) for i in range(N))
            mb.setObjective(totalb, GRB.MINIMIZE)
            
            B01a = ma.addMVar((M,M), vtype=GRB.BINARY, name="B01a")
            B10a = ma.addMVar((M,M), vtype=GRB.BINARY, name="B10a")
            B11a = ma.addMVar((M,M), vtype=GRB.BINARY, name="B11a")

            B01b = mb.addMVar((M,M), vtype=GRB.BINARY, name="B01b")
            B10b = mb.addMVar((M,M), vtype=GRB.BINARY, name="B10b")
            B11b = mb.addMVar((M,M), vtype=GRB.BINARY, name="B11b")

            z = mb.addMVar((1, N), vtype =GRB.BINARY, name="z")
            # Constraints (ensures no conflicts by checking each pair of columns (p, q)) 
            ma.addConstrs(-1 * X[i, p] + X[i, q] <= B01a[p,q] for p in range(M) for q in range(M) for i in range(N) if p != q)
            ma.addConstrs(X[i, p] - X[i, q] <= B10a[p,q] for p in range(M) for q in range(M) for i in range(N) if p != q)
            ma.addConstrs(X[i, p] + X[i, q] - 1 <= B11a[p,q] for p in range(M) for q in range(M) for i in range(N) if p != q)
            ma.addConstrs(B01a[p,q] + B10a[p,q] + B11a[p,q] <= 2 for p in range(M) for q in range(M) if p != q)

            mb.addConstrs(-1 * Y[i, p] + Y[i, q] <= B01b[p,q] for p in range(M) for q in range(M) for i in range(N) if p != q)
            mb.addConstrs(Y[i, p] - Y[i, q] <= B10b[p,q] for p in range(M) for q in range(M) for i in range(N) if p != q)
            mb.addConstrs(Y[i, p] + Y[i, q] - 1 <= B11b[p,q] for p in range(M) for q in range(M) for i in range(N) if p != q)
            mb.addConstrs(B01b[p,q] + B10b[p,q] + B11b[p,q] <= 2 for p in range(M) for q in range(M) if p != q)
        
            # Essential Partial Order Constraints
            mb.addConstr(z[i] <= (Y[u, i] - Y[v, i] + 1)/2 for i in range(N) for v in range(N))
            mb.addConstr(sum(z[i] for i in range(N)) >= 1)
            mb.optimize()
            
            if ma.objVal < mb.objVal:
                return True
            else: 
                return False

        except GurobiError as ex:
            print('*********ERROR*********')
            print(ex)

def Split(V):
    V = list(V)
    i = int(len(V)/2)
    return set(V[:i]), set(V[i:])

# Given a sample index u and set of test sample indices V, GetRelated(u, V) 
# outputs all samples v ∈ V such that u <e v
def GetRelated(u, V):
    if ILPincreased(u, V):
        if V.shape[0] == 1:
            return V
        else:
            Vl, Vr = Split(V)
            return GetRelated(u, Vl).union(GetRelated(u, Vr))
        
# Given a matrix D and positive integer k, GetEssential(D, k) outputs the 
# the <e relation
def GetEssential(D, k):
    S = {num+1 for num in range(D.shape[0])}
    ess_set = {}
    R = []
    for u in S:
        V = S.difference({u})
        R.append(GetRelated(u, V))
        ess_set = ess_set.union({(u, y) for y in R[u]})
    return ess_set


D = getMatrix("small_test.csv")
k = 10
print(GetEssential(D, k))


