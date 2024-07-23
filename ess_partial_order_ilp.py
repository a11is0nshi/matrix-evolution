"""
ILP Implementation inspired by (Malikic et. al, 867-868) paper. 
This implementation sets D to be a random binary m x n matrix, where there are 
m sequenced single cells and n mutations

sig1 is the minimum number of 0 to 1 flips necessary to make D a 
conflict free matrix while allowing for k 1 to 0 flips. 
"""
from gurobipy import *
from sys import * 
import numpy as np
import pandas as pd

# Change name and k to change file
name = "smallest_test.csv"
k = 1

D =  pd.read_csv(name).to_numpy()
n = D.shape[0]  # samples/rows
m = D.shape[1]  # mutations/cols

def addEssOrder(X, model, u, v, sigma):
    z = model.addMVar((m,), vtype=GRB.BINARY, name="z")
    model.addConstrs(X[u, i] - X[v, i] <= z[i] for i in range(m))
    model.addConstrs((X[u, i] - X[v, i] + 1) / 2 >= z[i] for i in range(m))
    model.addConstr(sum(z[i] for i in range(n)) >= 1)
    model.addConstr(sum(sum((1 - D[i, j])*(X[i, j]) for j in range(m)) for i in range(n)) <= sigma)

# Calculate sigma - minimum number of 0 -> 1 flips
def calcSigma(u, Vset, prime):
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

        total = sum(sum((1 - D[i, j])*(X[i, j]) for j in range(m)) for i in range(n))
        model.setObjective(total, GRB.MINIMIZE)

        # Enforce no conflicts by checking each pair of columns (p, q)) 
        model.addConstrs( X[i, q]- X[i, p] <= B01[p,q] for p in range(m) for q in range(p+1, m) for i in range(n))
        model.addConstrs(X[i, p] - X[i, q] <= B10[p,q] for p in range(m) for q in range(p+1, m) for i in range(n))
        model.addConstrs(X[i, p] + X[i, q] -1  <= B11[p,q] for p in range(m) for q in range(p+1, m) for i in range(n))
        model.addConstrs(B01[p,q] + B10[p,q] + B11[p,q] <= 2 for p in range(m) for q in range(p+1, m))
       
        # Ensure that there are at most k 1->0 flips
        model.addConstr(sum(D[i, j] * (1 - X[i, j]) for i in range(n) for j in range(m)) <= k)

        model.update()
        model.optimize()
        sigma = model.ObjVal

        if prime:
            for v in Vset:
                model = addEssOrder(X, model, u, v, sigma)
            model.update()
            model.optimize()
            return model.ObjVal
        else:
            return sigma

    

    except GurobiError as ex:
        print('*********ERROR*********')
        print(f"error: {ex}")

print(calcSigma(0, {1, 2}, False))

"""

def ILPincreased(u, Vset):
    sigma = calcSigma(u, Vset, False)
    sig_prime = calcSigma(u, Vset, True)
    print(f"u: {u}, V: {Vset}")
    print(f"    sigma: {sigma}, sig_prime: {sig_prime} \n")
    return sigma < sig_prime

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
        
# Given a matrix D and positive integer k, GetEssential(D, k) outputs <e
def GetEssential(D, k):
    S = {num for num in range(D.shape[0])}
    ess_set = set()
    R = []
    # for u in S:
    #     V = S.difference({u})
    #     R.append(GetRelated(u, V))
    #     P = {(u, y) for y in R[u]}
    #     temp = ess_set.union(P)
    #     ess_set = temp
    # return ess_set
    R.append(GetRelated(0, {1}))
    P = {(0, y) for y in R[0]}  
    return P
# print(f"R: {GetEssential(D, k)}")

"""