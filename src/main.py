# MAIN QMC PROGRAM:
import sys
import os
from datetime import datetime
from time import process_time

#-----------------------#
from energy import ENERGY
from qmc import QMC
from initialize import VARINIT
from pathdef import PATHS
from qgen import QGEN
from flibs import FORTLIBS
#-----------------------#
from openmm.app import *
from openmm import *
from openmm.unit import *

#-----------------------#
import numpy as np
import math
import pandas as pd
import random







#




#Begin Job:

job_start = datetime.now()
print("Begin Job at:",job_start)
jbegin = process_time()


#Instantiate class:
EN = ENERGY()
qmc = QMC()
vinit = VARINIT()
flib = FORTLIBS()
pths = PATHS()
qgen = QGEN()
#Get Input Vars:
vinit.read_inp('usermols.inp')
pths.fpaths()
#--------------------------------#
#Threshold:
thresh = 200.0
thresh2 = 1000.0


#Run Initial Sample:
if vinit.indq==1:
    s,b=qgen.sobol(vinit.npars,vinit.nstr)
qgen.writeang(s,vinit.npars,b)
#Initial Conformation generation:
flib.conf()

#Set files and folders:


#Use conformations and split them:
vinit.molsread(os.path.join(vinit.path,pths.molsfile),vinit.nstr,pths.spath)



#Loop over cycles to get structures:
#Applying Baker transformation:
mnum =[]
#open(fpdb,'w').close()
for i in (os.listdir(pths.spath)):
    pdb =   PDBFile(os.path.join(pths.spath,i))
    ENE = EN.qmcene(pdb)
    mno = int(i.partition('.')[0])


    
    print('****************************************************')
    print('MODEL NO:',mno)
    print("Intial Energy:",ENE)

    
    ncycles = 100
   
    if ENE>thresh2:
        print('RESAMPLING MOVES \n')
        bb=0
        while ENE>thresh2:
            bb = bb + 1
            print('CYCLE:',bb)
            ENE = qmc.RESAMPLE(ENE,mno)
            if bb>ncycles:
                break 
    
        print('----------------------------------------------------')
    if ENE>thresh and ENE<thresh2:
        print('APPLYING MOD M SHIFT TRANSFORMATION')
        cc=0
        while ENE>thresh:
            cc = cc + 1
            print('CYCLE:',cc)
            ENE = qmc.bakers_transform(ENE,mno)
            if cc > ncycles:
                break
        if ENE>thresh2:
            print('CONVERGENCE NOT MET: MOVING TO NEXT SAMPLE')         
            EN.outp(pdb,mno)
            break
    

        print('----------------------------------------------------')
    
    if ENE<thresh:
        
        qmc.accept_move(pdb,mno) 

       # print('****************************************************')
#Format conversions:
flib.pdbread()
#Cluster conformations:
flib.clust()
flib.dihs()



#qgen.writeang(s,npar,b)





job_end = datetime.now()
print("Job Ended at:",job_end)
jend = process_time()

print("Time Elapsed (s):",(jend-jbegin))
