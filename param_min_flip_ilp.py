"""
ILP Implementation inspired by (Malikic et. al, 867-868) paper. 
This implementation sets D to be a random binary m x n matrix, where there are 
m sequenced single cells and n mutations

"""
from gurobipy import *
from sys import * 
import numpy as np
import pandas as pd

def get_phyolin_matrix(name):
    df = pd.read_csv(name)
    return df.to_numpy()

try:
    # D is input binary matrix
    D = get_phyolin_matrix("Patient2_phyolin.csv")
    m, n = D.shape[0], D.shape[1]

    model = Model("min_flip_model")
  # X is a conflict free matrix by constraints
    X = model.addMVar((m,n), vtype=GRB.BINARY, name="X")

    # Objective function
    total = sum(sum(D[i, j]*(1 - X[i, j]) + (1 - D[i, j])*(X[i, j]) for j in range(n)) for i in range(m))
    model.setObjective(total, GRB.MINIMIZE)
    
    B01 = model.addMVar((n,n), vtype=GRB.BINARY, name="B01")
    B10 = model.addMVar((n,n), vtype=GRB.BINARY, name="B10")
    B11 = model.addMVar((n,n), vtype=GRB.BINARY, name="B11")

    z = model.addMVar((1, m), vtype =GRB.BINARY, name="z")
    u = model.addVar(vtype=GRB.INTEGER, name="u")
    v = model.addVar(vtype=GRB.INTEGER, name="v")
    # Constraints (ensures no conflicts by checking each pair of columns (p, q)) 
    model.addConstrs(-1 * X[i, p] + X[i, q] <= B01[p,q] for p in range(n) for q in range(n) for i in range(m) if p != q)
    model.addConstrs(X[i, p] - X[i, q] <= B10[p,q] for p in range(n) for q in range(n) for i in range(m) if p != q)
    model.addConstrs(X[i, p] + X[i, q] - 1 <= B11[p,q] for p in range(n) for q in range(n) for i in range(m) if p != q)
    model.addConstrs(B01[p,q] + B10[p,q] + B11[p,q] <= 2 for p in range(n) for q in range(n) if p != q)
   
   # Essential Partial Order Constraints
    model.addConstr(z[i] <= (X[u, i] - X[v, i] + 1)/2 for i in range(m))
    model.addConstr(sum(z[i] for i in range(m)) >= 1)
    model.addConstr(0 <= u < m)
    model.addConstr(0 <= v < m)
    model.addConstr(v != u)
    model.optimize()
    
    # Print results
    for y in model.getVars():
       print(y.varName, y.x)
    
    print(D)
    print('Optimal Objective function value:', model.objVal)

except GurobiError as ex:
    print('*********ERROR*********')
    print(ex)