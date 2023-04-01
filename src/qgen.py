import scipy
import os
import numpy as np
from scipy.stats import qmc
import math
import pandas as pd
import random
#Provide an index to choose Sobol-1/RandomLHS-2/OrhtogonalLHS-3/Halton-4 quasi-random sequences:


from initialize import VARINIT



vinit = VARINIT()
vinit.read_inp('usermols.inp')

class QGEN:

    def sobol(self,npar,nstr):
    #For sobol instr to be converted to base 2:
        b = int(math.log(nstr,2))    
        so = qmc.Sobol(npar)
        s_so=so.random_base2(m=b)
        s_sa = s_so*vinit.search_sp
        sample = s_sa
        print("Randomised Sobol sequence generated\n")
        print("Discrepancy:",qmc.discrepancy(s_so))
   
        return sample,2**b
    


    def rshift(self,sample,nstr,npar,search_sp):
        #Shift nsamples to base 2:
        b = int(math.log(nstr,2))
        #Choose a uniform random variable:
        U=np.zeros((nstr,npar))
        for i in range(nstr):
            U[i] = random.uniform(0.0,search_sp)
        
        U = np.around(U,2) 
        #Use shift:
        s_shifted = np.zeros((nstr,npar))
        for i in range(nstr):
            for j in range(npar):
                s_shifted[i,j]=sample[i,j]-U[i,j]
    
        return s_shifted



#def restep(search_sp)


    def ord_sobol(self,npar,nstr):
        #Convert to base 2:
        b = int(math.log(nstr,2))
        so = qmc.Sobol(npar,scramble=False)
        s_so=so.random_base2(m=b)
        s_sa = s_so*vinit.search_sp
        sample = s_sa
        print("Ordered Sobol sequence generated\n")
        print("Discrepancy:",qmc.discrepancy(s_so))
   
        return sample,2**b
  




    def writeang(self,sz,npar,b):
        # en = np.zeros(b)
        # sz = np.column_stack((en,sample))
        sz = np.around(sz,decimals=2)
        df = pd.DataFrame(sz)
        df.to_csv('mols_q.out',header=None,sep='\t',index=True)




#search_sp = float(input('Enter the search space in degrees(Ex:180.0):'))


#nstr,ind,npar,path,mname,search_sp  =  read_inp('usermols.inp')
#if ind==1:
#    s,b=sobol(npar,nstr)
#if ind==2:
#    s,b=ord_sobol(npar,nstr)
#writeang(s,npar,b)
#l = rshift(s,nstr,npar,search_sp)
#print(s)
#print("Shifted")
#print(l)
#if ind==3:
#    s,b=mols(npar,nstr)
#if ind==4:
#   s,b=ord_sobol(npar,nstr)
#if ind==5:
#    s,b = uniform(npar,nstr)
#writeang(s,npar,b)



