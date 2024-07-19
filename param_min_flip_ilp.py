"""
ILP Implementation inspired by (Malikic et. al, 867-868) paper. 
This implementation sets D to be a random binary m x n matrix, where there are 
m sequenced single cells and n mutations

"""
from gurobipy import *
from sys import * 
import numpy as np
import pandas as pd



def getMatrix(name):
    df = pd.read_csv(name)
    return df.to_numpy()

# given a Matrix M and an index u, getV(M, u) returns V = M \ M[u]
def getV(M, u):
    return np.delete(M, u, 0)


def ILPincreased(u, V):

    try:
        # D is input binary matrix
        # D = get_phyolin_matrix("Patient2_phyolin.csv")
        m, n = V.shape[0], V.shape[1]

        m1 = Model("min_flip_model1")
        m2 = Model("min_flip_model2")

        
    # X is a conflict free matrix by constraints
        X = m1.addMVar((m,n), vtype=GRB.BINARY, name="X")
        Y = m2.addMVar((m,n), vtype=GRB.BINARY, name="Y")




        # Objective function
        total = sum(sum(V[i, j]*(1 - X[i, j]) + (1 - V[i, j])*(X[i, j]) for j in range(n)) for i in range(m))
        m1.setObjective(total, GRB.MINIMIZE)

        total2 = sum(sum(V[i, j]*(1 - Y[i, j]) + (1 - V[i, j])*(Y[i, j]) for j in range(n)) for i in range(m))
        m2.setObjective(total, GRB.MINIMIZE)
        
        B01 = m1.addMVar((n,n), vtype=GRB.BINARY, name="B01")
        B10 = m1.addMVar((n,n), vtype=GRB.BINARY, name="B10")
        B11 = m1.addMVar((n,n), vtype=GRB.BINARY, name="B11")

        B012 = m1.addMVar((n,n), vtype=GRB.BINARY, name="B01")
        B102 = m1.addMVar((n,n), vtype=GRB.BINARY, name="B10")
        B112 = m1.addMVar((n,n), vtype=GRB.BINARY, name="B11")

        z = m2.addMVar((1, m), vtype =GRB.BINARY, name="z")
        u = m2.addVar(vtype=GRB.INTEGER, name="u")
        v = m2.addVar(vtype=GRB.INTEGER, name="v")
        V = m2.addVar(vtype=GRB.INTEGER, name="V")

        # Constraints (ensures no conflicts by checking each pair of columns (p, q)) 
        m1.addConstrs(-1 * X[i, p] + X[i, q] <= B01[p,q] for p in range(n) for q in range(n) for i in range(m) if p != q)
        m1.addConstrs(X[i, p] - X[i, q] <= B10[p,q] for p in range(n) for q in range(n) for i in range(m) if p != q)
        m1.addConstrs(X[i, p] + X[i, q] - 1 <= B11[p,q] for p in range(n) for q in range(n) for i in range(m) if p != q)
        m1.addConstrs(B01[p,q] + B10[p,q] + B11[p,q] <= 2 for p in range(n) for q in range(n) if p != q)

        m2.addConstrs(-1 * Y[i, p] + Y[i, q] <= B012[p,q] for p in range(n) for q in range(n) for i in range(m) if p != q)
        m2.addConstrs(Y[i, p] - Y[i, q] <= B102[p,q] for p in range(n) for q in range(n) for i in range(m) if p != q)
        m2.addConstrs(Y[i, p] + Y[i, q] - 1 <= B112[p,q] for p in range(n) for q in range(n) for i in range(m) if p != q)
        m2.addConstrs(B012[p,q] + B102[p,q] + B112[p,q] <= 2 for p in range(n) for q in range(n) if p != q)
    
    # Essential Partial Order Constraints
        m2.addConstr(z[i] <= (Y[u, i] - Y[v, i] + 1)/2 for i in range(m))
        m2.addConstr(sum(z[i] for i in range(m)) >= 1)
        m2.addConstr(0 <= u)
        m2.addConstr(u < m)
        m2.addConstr(0 <= v)
        m2.addConstr(v < m)
        m2.addConstr(v != u)
        m2.optimize()
        
        # # Print results
        # for y in m1.getVars():
        # print(y.varName, y.x)
        
        if m1.objVal < m2.objVal:
            return True
        else: 
            return False

    except GurobiError as ex:
        print('*********ERROR*********')
        print(ex)

def Split(V):
    center = V.shape[0] // 2
    Vl = np.zeros(len(center), V.shape[1])
    Vr = np.zeros(V.shape[0]-len(center), V.shape[1])
    for i in range(center):
        Vl.append(V[i])
    for j in range(V.shape[0]-len(center)):
        Vr.append(V[j])
    return V, V

def GetRelated(u, V):
    if ILPincreased:
        if V.shape[0] == 1:
            return V
        else:
            (Vl, Vr) = Split(V)
            return np.concatenate(GetRelated(u, Vl), GetRelated(u, Vr))
        

M = getMatrix("nameoffile")
Ru = np.empty(M.shape[0], dtype=object)
for u in range(M.shape[0]):
    V = getV(M, u)
    Ru[u] = GetRelated(u, V)
print(str(Ru))

