"""
ILP Implementation inspired by (Malikic et. al, 867-868) paper. 
This implementation sets D to be a random binary m x n matrix, where there are 
m sequenced single cells and n mutations

"""
from gurobipy import *
import numpy as np
import process_data as p
import pandas as pd

# Change name and k to change file
name = "smallest_test.csv"

D =  pd.read_csv(name).to_numpy()

try:
    # D is input binary matrix
    m, n = D.shape[0], D.shape[1]

    model = Model("min_flip_model")

    # X is a conflict free matrix by constraints
    X = model.addMVar((m,n), vtype=GRB.BINARY, name="X")

    # Objective function
    total = sum(sum((1 - D[i, j])*(X[i, j]) for j in range(n)) for i in range(m))
    model.setObjective(total, GRB.MINIMIZE)

    k = 0
    model.addConstr(sum(D[i, j] * (1 - X[i, j]) for i in range(m) for j in range(n)) <= k)
    
    B01 = model.addMVar((n,n), vtype=GRB.BINARY, name="B01")
    B10 = model.addMVar((n,n), vtype=GRB.BINARY, name="B10")
    B11 = model.addMVar((n,n), vtype=GRB.BINARY, name="B11")

    # Constraints (ensures no conflicts by checking each pair of columns (p, q)) 
    model.addConstrs(-1 * X[i, p] + X[i, q] <= B01[p,q] for p in range(n) for q in range(n) for i in range(m) if p != q)
    model.addConstrs(X[i, p] - X[i, q] <= B10[p,q] for p in range(n) for q in range(n) for i in range(m) if p != q)
    model.addConstrs(X[i, p] + X[i, q] - 1 <= B11[p,q] for p in range(n) for q in range(n) for i in range(m) if p != q)
    model.addConstrs(B01[p,q] + B10[p,q] + B11[p,q] <= 2 for p in range(n) for q in range(n) if p != q)
    
    model.optimize()
    print(f"sigma: {model.ObjVal}")
    
    # f = open("Patient_6_Results.txt", "w")
    # # Print results
    # f.write("D")
    # D_str = np.array2string(D, threshold = np.inf)
    # f.write(f"\n{D_str} \n \n")
    # f.write("X")
    # values = model.getAttr("X", model.getVars())
    # values = np.array(values).reshape((m, n))
    # X_str = np.array2string(X.X, threshold = np.inf)
    # f.write(f"\n{X_str} \n \n")
    # f.write("Optimal Objective function value: \n")
    # f.write(str(model.objVal))
    # f.close()


except GurobiError as ex:
    print('*********ERROR*********')
    print(ex)