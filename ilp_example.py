"""
Attempt to implement ILP example from Dr. Yaw's slides
"""

from gurobipy import *

try:
   # Create a new model
   m = Model("mip1")
   # Create variables
   x = m.addVar(vtype=GRB.BINARY, name="x")
   y = m.addVar(vtype=GRB.BINARY, name="y")
    
   # Set objective
   m.setObjective(5 * x + 8 * y, GRB.MAXIMIZE)

   # Add constraint:
   m.addConstr(x + y <= 6, "c0")
   # Add constraint:
   m.addConstr(5 * x + 9 * y <= 45, "c1")
   # Add constraint:
   m.addConstr(x >= 0, "c2")
   m.addConstr(y >= 0, "c3")

   # m.update()
   m.optimize()

   for v in m.getVars():
      print(v.varName, v.x)
   
   print('Optimal Objective function value:', m.objVal)

except GurobiError:
    print('Error reported')