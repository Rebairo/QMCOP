#Convert the points distributed using sobol sequence in the unit hypercube to angles.
import numpy as np
import os
import pandas as pd
import math
import qgen

d = pd.read_csv("mols_q.out",header=None,delim_whitespace=True)
ds = d.iloc[:,2]

nstr,ind,npar = qgen.read_inp("usermols.inp")
bits = int(math.log(nstr,2))

def dtb(n):
    return bin(n).replace("0b", "")

x=[]

for i,j in enumerate(ds):
    x.append(j*360.0)
print()
