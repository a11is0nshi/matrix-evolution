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

def ILP(u, Vset, prime):
    N = D.shape[0]  # num of samples - rows
    M = D.shape[1]  # num of mutations - cols

    try:
        m = Model("min_flip_model")
    
        # X is a conflict free matrix by constraints
        X = m.addMVar((N, M), vtype=GRB.BINARY, name="X")

        # Objective function
        total = sum(sum((1 - D[i, j])*(X[i, j]) for j in range(M)) for i in range(N))
        m.setObjective(total, GRB.MINIMIZE)
        
        # Add auxillary constraints
        B01 = m.addMVar((M,M), vtype=GRB.BINARY, name="B01")
        B10 = m.addMVar((M,M), vtype=GRB.BINARY, name="B10")
        B11 = m.addMVar((M,M), vtype=GRB.BINARY, name="B11")

        # Constraints (ensures no conflicts by checking each pair of columns (p, q)) 
        m.addConstrs(-1 * X[i, p] + X[i, q] <= B01[p,q] for p in range(M) for q in range(M) for i in range(N) if p != q)
        m.addConstrs(X[i, p] - X[i, q] <= B10[p,q] for p in range(M) for q in range(M) for i in range(N) if p != q)
        m.addConstrs(X[i, p] + X[i, q] - 1 <= B11[p,q] for p in range(M) for q in range(M) for i in range(N) if p != q)
        m.addConstrs(B01[p,q] + B10[p,q] + B11[p,q] <= 2 for p in range(M) for q in range(M) if p != q)

        # Essential Partial Order Constraints
        if prime:
            z = m.addMVar((N,), vtype =GRB.BINARY, name="z")
            print(z)
            m.addConstrs(z[i] <= (X[u-1, i] - X[v-1, i] + 1)/2 for i in range(N) for v in Vset)
            m.addConstr(sum(z[i] for i in range(N)) >= 1)

        m.optimize()
        return m.ObjVal

    except GurobiError as ex:
        print('*********ERROR*********')
        print(ex)

def ILPincreased(u, Vset):
    r1 = ILP(u, Vset, False)
    r2 = ILP(u, Vset, True)
    return r1 < r2

def Split(V):
    V = list(V)
    i = int(len(V)/2)
    return set(V[:i]), set(V[i:])

# Given a sample index u and set of test sample indices V, GetRelated(u, V) 
# outputs all samples v ∈ V such that u <e v
def GetRelated(u, V):
    if ILPincreased(u, V):
        if len(V) == 1:
            return V
        else:
            Vl, Vr = Split(V)
            return GetRelated(u, Vl).union(GetRelated(u, Vr))
        
# Given a matrix D and positive integer k, GetEssential(D, k) outputs the 
# the <e relation
def GetEssential(D, k):
    S = {num+1 for num in range(D.shape[0])}
    print(f"S : {S}")
    ess_set = set()
    R = []
    for u in S:
        V = S.difference({u})
        R.append(GetRelated(u, V))
        if R[-1] == None:
            print(f"u: {u}, V: {V}")
        print(F"R : {R}")
        P = {(u, y) for y in R[u-1]}
        temp = ess_set.union(P)
        ess_set = temp
    return ess_set

D = getMatrix("small_test.csv")
k = 10
print(f"GetEssential: {GetEssential(D, k)}")


