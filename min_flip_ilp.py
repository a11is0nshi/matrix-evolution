"""
ILP Implementation inspired by (Malikic et. al, 867-868) paper
"""
from gurobipy import *
import numpy as np

try:
    m = 6
    n = 4
    # random binary matrix of size m x n
    # rows correspond to sequenced single cells, columns to mutations
    D = np.random.choice(2, m*n).reshape((m, n))

    # Create a new model
    model = Model("min_flip_model")

    # X is constrained to be a conflict free matrix
    X = model.addMVar((m,n), vtype=GRB.BINARY, name="X")

    # Set objective function
    total = sum(sum(D[i, j]*(1 - X[i, j]) + (1 - D[i, j])*(X[i, j]) for j in range(n)) for i in range(m))
    model.setObjective(total, GRB.MINIMIZE)
    
    """
    # Initialize values of X
    for i in range(m):
        for j in range(n):
            # Create variables (entries of X)
            model.addVar(vtype=GRB.BINARY, name=f"X[{i}, {j}]")
    """
    
    B01 = model.addMVar((n,n), vtype=GRB.BINARY, name="B01")
    B10 = model.addMVar((n,n), vtype=GRB.BINARY, name="B10")
    B11 = model.addMVar((n,n), vtype=GRB.BINARY, name="B11")

    model.update()

    # Add constraints to ensure that there is no conflict involving cols p and q. 
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