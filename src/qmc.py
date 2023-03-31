# Compile the files using f2py and import as libraries here:



import numpy as np
import pandas as pd
import os
import energy
import qgen
import initialize
import pathdef
import flibs
from openmm.app import *
from openmm import *
from openmm.unit import *


#Initialize VARS:

nstr,ind,npar,path,mname,search_sp  =  initialize.read_inp('usermols.inp')

spath,mpath,molsfile = pathdef.path(path,mname)
mnum=[]
def RESAMPLE(ENE,mno):
    
    mnum.append(mno)
        
    with open ('mols_q.out','r+') as ff:
        for j in ff.readlines():
            if(mno==int(j.split()[0])):
                hh = j
                #print(hh) 
            else:
                hh = j
    jj = hh.split()[1:]
        #udat = np.around(np.asarray(jj,dtype=float),decimals=2)
    U2,b = qgen.sobol(npar,1) 
        #sdat =(((udat-(1-2*U))%2))
    s2dat = np.transpose(U2)
    s2dat2 = s2dat.reshape((1,npar))
    df2_udat = pd.DataFrame(s2dat2)
    df2_udat.to_csv("mols_q.out",sep="\t",header=None,index=True)
    flibs.conf()
        
    pdb =   PDBFile(os.path.join(path,molsfile))
    ENE = energy.energy(pdb)
    print('ENERGY AFTER RESAMPLING:',ENE)

    return ENE

def bakers_transform(ENE,mno):

    mnum.append(mno)
        
    with open ('mols_q.out','r+') as ff:
        for j in ff.readlines():
            if(mno==int(j.split()[0])):
                hh = j
            else:
                hh = j
    jj = hh.split()[1:]
    udat = np.around(np.asarray(jj,dtype=float),decimals=2)
    U = (np.random.random(npar))*(search_sp)
    sdat =(((udat-(U))))
    sdat2 = sdat.reshape((1,npar))
    df_udat = pd.DataFrame(sdat2)
    df_udat.to_csv("mols_q.out",sep="\t",header=None,index=True)
    flibs.conf()
    
    pdb =   PDBFile(os.path.join(path,molsfile))
    ENE = energy.energy(pdb)
    print('ENERGY AFTER MOD SHIFT TRANSFORMATION:',ENE)
            
    return ENE


def accept_move(pdb,mname,mno,mpath,path):
    
    print('MINIMISING OPTIMAL CONFORMATION')
    energy.outp(pdb,mname,mno,mpath,path)
    print('OPTIMAL CONFORMATION GENERATED')
