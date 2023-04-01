
# Compile the files using f2py and import as libraries here:



import numpy as np
import pandas as pd
import os
from energy import ENERGY
from qgen import QGEN
from initialize import VARINIT
from pathdef import PATHS
from flibs import FORTLIBS
from openmm.app import *
from openmm import *
from openmm.unit import *


#Initialize VARS:

EN = ENERGY()
vinit = VARINIT()
flib = FORTLIBS()
pths = PATHS()
qgen = QGEN()
inpf = vinit.read_inp('usermols.inp')
pths.fpaths()
mnum=[]


class QMC:
    def RESAMPLE(self,ENE,mno):
    
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
        U2,b = qgen.sobol(vinit.npars,1) 
        #sdat =(((udat-(1-2*U))%2))
        s2dat = np.transpose(U2)
        s2dat2 = s2dat.reshape((1,vinit.npars))
        df2_udat = pd.DataFrame(s2dat2)
        df2_udat.to_csv("mols_q.out",sep="\t",header=None,index=True)
        flib.conf()
        
        pdb =   PDBFile(os.path.join(vinit.path,pths.molsfile))
        ENE = EN.qmcene(pdb)
        print('ENERGY AFTER RESAMPLING:',ENE)

        return ENE

    def bakers_transform(self,ENE,mno):

        mnum.append(mno)
        
        with open ('mols_q.out','r+') as ff:
            for j in ff.readlines():
                if(mno==int(j.split()[0])):
                    hh = j
                else:
                    hh = j
        jj = hh.split()[1:]
        udat = np.around(np.asarray(jj,dtype=float),decimals=2)
        U = (np.random.random(vinit.npars))*(vinit.search_sp)
        sdat =(((udat-(U))))
        sdat2 = sdat.reshape((1,vinit.npars))
        df_udat = pd.DataFrame(sdat2)
        df_udat.to_csv("mols_q.out",sep="\t",header=None,index=True)
        flib.conf()
    
        pdb =   PDBFile(os.path.join(vinit.path,pths.molsfile))
        ENE = EN.qmcene(pdb)
        print('ENERGY AFTER MOD SHIFT TRANSFORMATION:',ENE)
            
        return ENE


    def accept_move(self,pdb,mno):
    
        print('MINIMISING OPTIMAL CONFORMATION')
        EN.outp(pdb,vinit.mname,mno,pths.mpath,vinit.path)
        print('OPTIMAL CONFORMATION GENERATED')
