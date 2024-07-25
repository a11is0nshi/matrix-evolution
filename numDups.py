import pandas as pd
import numpy as np

name = "AML-10_rep10.csv"

df = pd.read_csv(name)
D = df.to_numpy() # has duplicates
df.drop_duplicates(inplace=True)
E = df.to_numpy(dtype=int) # no duplicates

N = {}
# print(D)
# print(E)

for i in range(E.shape[0]):
    dups = 0
    for x in range(D.shape[0]):
        if np.all(D[x] == E[i]):
            dups += 1
    N[tuple(E[i])] = dups
print(N)
