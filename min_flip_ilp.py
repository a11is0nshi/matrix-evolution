"""
ILP Implementation inspired by (Malikic et. al, 867-868) paper. 
This implementation sets D to be a random binary m by n matrix, where there are 
m sequenced single cells and n mutations

"""
from gurobipy import *
import numpy as np

try:
    # D random binary matrix of size m x n
    m = 6
    n = 4
    D = np.random.choice(2, m*n).reshape((m, n))

    model = Model("min_flip_model")

    # X is coflict free matrix by constraints
    X = model.addMVar((m,n), vtype=GRB.BINARY, name="X")

    # Objective function
    total = sum(sum(D[i, j]*(1 - X[i, j]) + (1 - D[i, j])*(X[i, j]) for j in range(n)) for i in range(m))
    model.setObjective(total, GRB.MINIMIZE)
    
    B01 = model.addMVar((n,n), vtype=GRB.BINARY, name="B01")
    B10 = model.addMVar((n,n), vtype=GRB.BINARY, name="B10")
    B11 = model.addMVar((n,n), vtype=GRB.BINARY, name="B11")

    model.update()

    # Constraints (ensures no conflicts by checking each pair of columns (p, q)) 
    model.addConstrs(-1 * X[i, p] + X[i, q] <= B01[p,q] for p in range(n) for q in range(n) for i in range(m) if p != q)
    model.addConstrs(X[i, p] - X[i, q] <= B10[p,q] for p in range(n) for q in range(n) for i in range(m) if p != q)
    model.addConstrs(X[i, p] + X[i, q] - 1 <= B11[p,q] for p in range(n) for q in range(n) for i in range(m) if p != q)
    model.addConstrs(B01[p,q] + B10[p,q] + B11[p,q] <= 2 for p in range(n) for q in range(n) if p != q)
    
    model.optimize()
    
    # Print results
    for v in model.getVars():
       print(v.varName, v.x)
    
    print(D)
    print('Optimal Objective function value:', model.objVal)

except GurobiError as ex:
    print('*********ERROR*********')
    print(ex)