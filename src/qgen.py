import scipy
import os
import numpy as np
from scipy.stats import qmc
import math
import pandas as pd
import random
#Provide an index to choose Sobol-1/RandomLHS-2/OrhtogonalLHS-3/Halton-4 quasi-random sequences:

def read_inp(inpfile):

    with open(inpfile,'r')as f:
        inp=[]
        for i,j in enumerate(f):
            if i==0:
                path = str(j.split()[0])
            if i==1:
                mname = str(j.split()[0])
            if i==3:
                nstr=int(j.split()[0])
            if i==10:
                indq=int(j.split()[1])
            if i==11:
                npars=int(j.split()[1])
            if i==12:
                search_sp=float(j.split()[1])
          
                
    return nstr,indq,npars,path,mname,search_sp



def sobol(npar,nstr):
    #For sobol instr to be converted to base 2:
    b = int(math.log(nstr,2))    
    so = qmc.Sobol(npar)
    s_so=so.random_base2(m=b)
    s_sa = s_so*search_sp
    sample = s_sa
    print("Randomised Sobol sequence generated\n")
    print("Discrepancy:",qmc.discrepancy(s_so))
   
    return sample,2**b
    


def rshift(sample,nstr,npar,search_sp):
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


def ord_sobol(npar,nstr):
    #Convert to base 2:
    b = int(math.log(nstr,2))
    so = qmc.Sobol(npar,scramble=False)
    s_so=so.random_base2(m=b)
    s_sa = s_so*search_sp
    sample = s_sa
    print("Ordered Sobol sequence generated\n")
    print("Discrepancy:",qmc.discrepancy(s_so))
   
    return sample,2**b
  



# def rlhs(npar,nstr):
    # rlhs=qmc.LatinHypercube(d=npar,optimization='random-cd')
    # r_s = rlhs.random(n=nstr)
    # r_s = np.around(r_s,decimals=2)
    # rs = r_s*search_sp
    # sample = rs
    # print("Randomised Latin Squares sequences generated \n")
    # print("Discrepancy:",qmc.discrepancy(r_s))
    # return sample,nstr
  

# def mols(npar,nstr):
    # if npar>math.sqrt(nstr):
        # print("Error: Ntors should be lesser than sqrt of no_of_str")
        # exit()
    # mls=qmc.LatinHypercube(d=npar,strength=2)
    # k_sq=math.sqrt(nstr)
    # #for i in range(2,nstr):
    # #    if (k_sq%i)==0:
    # #        print("Error: Enter No of structures = square of a prime number ")
    # #        exit()
        
    # mo_ls=mls.random(n=nstr)
    # m_ls = np.around(mo_ls,decimals=2)
    # ms = m_ls*search_sp
    # sample = ms
    # print("Orthogonal Latin square sequences generated\n")
    # print("Discrepancy:",qmc.discrepancy(m_ls))
    # return sample,nstr


# def uniform(npar,nstr):
    # s = np.zeros((nstr,npar))
    # for i in range(nstr):
# def uniform(npar,nstr):
    # s = np.zeros((nstr,npar))
    # for i in range(nstr):
 
        # s[i]=(np.random.uniform(npar))
    # sample = s*search_sp
    # return sample,nstr



def writeang(sz,npar,b):
    # en = np.zeros(b)
    # sz = np.column_stack((en,sample))
     sz = np.around(sz,decimals=2)
     df = pd.DataFrame(sz)
     df.to_csv('mols_q.out',header=None,sep='\t',index=True)




#search_sp = float(input('Enter the search space in degrees(Ex:180.0):'))


nstr,ind,npar,path,mname,search_sp  =  read_inp('usermols.inp')
if ind==1:
    s,b=sobol(npar,nstr)
if ind==2:
    s,b=ord_sobol(npar,nstr)
writeang(s,npar,b)
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



