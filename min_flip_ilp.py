"""
Implementation of ILP from Tumor Evolution paper (Malikic et. al, 867-868)
This ILP formulation searches for the conflict-free matrix X for which 
P(D | X) is maximized. 
"""

from gurobipy import *
import numpy as np

# return 1 if there exists a row i such that X[i, p] = a, X[i, q] = b
def b_func(X, p, q, a, b):
    # extract number of rows in matrix X
    m = X.shape[0]
    for i in range(m):
        if X[i, p] == a and X[i, q] == b:
            return 1
    return 0

try:
    m = 6
    n = 4
    # m x n binary matrix from random 1-d array of size 24 from np.arange(2)
    # rows correspond to sequenced single cells, columns to mutations
    D = np.random.choice(2, 24).reshape((m, n))

    # X is constrained to be a conflict free matrix
    X = np.zeros((m, n), dtype=object)

    # Create a new model
    model = Model("mip1")
    
    # Initialize values of X
    for i in range(m):
        for j in range(n):
            # Create variables (entries of X)
            x_i_j = model.addVar(vtype=GRB.BINARY, name=f"X[{i}, {j}]")
            X[i, j] = x_i_j

    # Add constraints to ensure that there is no conflict involving cols p and q. 
    count = 0
    for p in range(n):
        for q in range(p+1, n):
            for i in range(m):
                model.addConstr(b_func(X, p, q, 0, 1) >= (-1 * X[i, p] + X[i, q]), f"c{count}")
                count += 1
                model.addConstr(b_func(X, p, q, 1, 0) >= (X[i, p] - X[i, q]), f"c{count}")
                count += 1
                model.addConstr(b_func(X, p, q, 1, 1) >= (X[i, p] + X[i, q] - 1), f"c{count}")
                count += 1
                model.addConstr(b_func(X, p, q, 0, 1) + b_func(X, p, q, 1, 0) + b_func(X, p, q, 1, 1) <= 2, f"c{count}")
                count += 1
    
    # Set objective
    total = 0
    for i in range(m):
        for j in range(n):
            total += D[i, j]*(1 - X[i, j]) + (1 - D[i, j])*(1 - X[i, j]) + (1 - D[i, j])*(X[i, j]) + D[i, j]*(X[i, j])
    
    model.setObjective(total, GRB.MAXIMIZE)

    # Print results
    for v in model.getVars():
        print(v.varName, v.x)
   
    print('Optimal Objective function value:', model.objVal)

except GurobiError as ex:
    print('ERROR')
    print(ex)