"""
ILP Implementation inspired by (Malikic et. al, 867-868) paper
"""
from gurobipy import *
import numpy as np

# # returns 1 if there exists a row i such that X[i, p] = a, X[i, q] = b
# def b_func(X, p, q, a, b):
#     # extract number of rows in matrix X
#     m = X.shape[0]
#     for i in range(m):
#         # !!! X contains gurobipy.Var while a, b are ints !!!
#         if X[i, p] == a and X[i, q] == b:
#             return 1
#     return 0

try:
    m = 6
    n = 4
    # random binary matrix of size m x n
    # rows correspond to sequenced single cells, columns to mutations
    D = np.random.choice(2, m*n).reshape((m, n))

    # Create a new model
    model = Model("mip1")

    # X is constrained to be a conflict free matrix
    X = model.addMVar((m,n), vtype=GRB.BINARY)
    
    # Initialize values of X
    for i in range(m):
        for j in range(n):
            # Create variables (entries of X)
            model.addVar(vtype=GRB.BINARY, name=f"X[{i}, {j}]")
    
    B01 = np.zeros((n, n), dtype=object)
    B10 = np.zeros((n, n), dtype=object)
    B11 = np.zeros((n, n), dtype=object)
    for p in range(n):
        for q in range(p+1, n):
            B01[p, q] = model.addVar(vtype=GRB.BINARY, name=f"B[{p}, {q}, 0, 1]")
            B10[p, q] = model.addVar(vtype=GRB.BINARY, name=f"B[{p}, {q}, 1, 0]")
            B11[p, q] = model.addVar(vtype=GRB.BINARY, name=f"B[{p}, {q}, 1, 1]")

    # Add constraints to ensure that there is no conflict involving cols p and q. 
    count = 0
    for p in range(n):
        for q in range(p+1, n):
            for i in range(m):
                model.addConstr(-X[i, p] + X[i, q] <= B01[p,q], f"c{count}")
                count += 1
                model.addConstr(X[i, p] - X[i, q] <= B10[p,q], f"c{count}")
                count += 1
                model.addConstr(X[i, p] + X[i, q] - 1 <= B11[p,q], f"c{count}")
                count += 1
                model.addConstr(B01[p,q] + B10[p,q] + B11[p,q] <= 2, f"c{count}")
                count += 1
    
    # Set objective function - equal to num of bit flips?    
   
    total = quicksum(quicksum(D[i, j]*(1 - X[i, j]) + (1 - D[i, j])*(1 - X[i, j]) + (1 - D[i, j])*(X[i, j]) + D[i, j]*(X[i, j]) for j in range(n)) for i in range(m))
    
    model.setObjective(total, GRB.MINIMIZE)
    model.optimize()
    # Print results
    #for v in model.getVars():
       # print(v.varName, v.x)
    
    print(D)
    print('Optimal Objective function value:', model.getObjective())

except GurobiError as ex:
    print('*********ERROR*********')
    print(ex)