"""
ILP Implementation inspired by (Malikic et. al, 867-868) paper. 
This implementation sets D to be a random binary m x n matrix, where there are 
m sequenced single cells and n mutations
"""
from gurobipy import *
from sys import * 
import numpy as np
import pandas as pd

# Change name and k to change file
name = "small_test.csv"
k = 10

D =  pd.read_csv(name).to_numpy()

def ILP(u, Vset, prime):
    N = D.shape[0]  # num of samples - rows
    M = D.shape[1]  # num of mutations - cols
    Nv = len(Vset)

    try:
        env = Env(empty=True)
        env.setParam("OutputFlag",0)
        env.start()
        m = Model("min_flip_model", env = env)
        m.Params.LogToConsole = 0
    
        # X is a conflict free matrix by constraints
        X = m.addMVar((N, M), vtype=GRB.BINARY, name="X")

        # Objective function
        total = sum(sum((1 - D[i, j])*(X[i, j]) for j in range(M)) for i in range(N))
        m.setObjective(total, GRB.MINIMIZE)
        
        # Add auxillary constraints
        B01 = m.addMVar((M, M), vtype=GRB.BINARY, name="B01")
        B10 = m.addMVar((M, M), vtype=GRB.BINARY, name="B10")
        B11 = m.addMVar((M, M), vtype=GRB.BINARY, name="B11")

        # Constraints (ensures no conflicts by checking each pair of columns (p, q)) 
        m.addConstrs(X[i, q] - X[i, p] <= B01[p,q] for p in range(M) for q in range(p+1, M) for i in range(N))
        m.addConstrs(X[i, p] - X[i, q] <= B10[p,q] for p in range(M) for q in range(p+1, M) for i in range(N))
        m.addConstrs(X[i, p] + X[i, q] <= B11[p,q] for p in range(M) for q in range(p+1, M) for i in range(N))
        m.addConstrs(B01[p,q] + B10[p,q] + B11[p,q] <= 2 for p in range(M) for q in range(p+1, M))
        m.addConstr(sum(D[i, j] * (1 - X[i, j]) for i in range(N) for j in range(M)) <= k)
        m.update()
        m.optimize()

        # Essential Partial Order Constraints
        if prime:
            z = m.addMVar((M,Nv), vtype=GRB.BINARY, name="z")
            m.update()
            for i in range(M):
                for v in range(Nv):
                    m.addConstr(z[i, v] <= (X[u-1, i] - X[list(Vset)[v]-1, i] + 1)/2)
            
            m.addConstr(sum(z[i, v] for i in range(M) for v in range(Nv)) <= 1)
            print(f"prime: {m.NumConstrs} constrs and {m.NumBinVars} vars")

        else:
            print(f"initial: {m.NumConstrs} constrs and {m.NumBinVars} vars")

        m.optimize()
        m.write('model.mps')
        toReturn = m.ObjVal
        m.reset()
        return toReturn

    except GurobiError as ex:
        print('*********ERROR*********')
        print(f"error: {ex}")

def ILPincreased(u, Vset):
    sig1 = ILP(u, Vset, False)
    sig2 = ILP(u, Vset, True)
    print(f"u: {u}, V: {Vset}")
    print(f"    sig1: {sig1}, sig2: {sig2} \n")
    return sig1 < sig2

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
    else: 
        return {}
        
# Given a matrix D and positive integer k, GetEssential(D, k) outputs the 
# the <e relation
def GetEssential(D, k):
    S = {num+1 for num in range(D.shape[0])}
    ess_set = set()
    R = []
    for u in S:
        V = S.difference({u})
        R.append(GetRelated(u, V))
        P = {(u, y) for y in R[u-1]}
        temp = ess_set.union(P)
        ess_set = temp
    return ess_set


print(f"R: {GetEssential(D, k)}")


