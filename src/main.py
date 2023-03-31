# MAIN QMC PROGRAM:
import sys
import os
from datetime import datetime
from time import process_time

#-----------------------#
import energy
import qmc
import initialize
import pathdef
import qgen
import flibs
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


#Get Input Vars:
nstr,ind,npar,path,mname,search_sp  =  initialize.read_inp('usermols.inp')

#--------------------------------#
#Threshold:
thresh = 200.0
thresh2 = 1000.0


#Run Initial Sample:
if ind==1:
    s,b=qgen.sobol(npar,nstr)

#Initial Conformation generation:
flibs.conf()

#Set files and folders:

spath,mpath,molsfile = pathdef.path(path,mname)

#Use conformations and split them:
sp = initialize.molsread(os.path.join(path,molsfile),nstr,spath)



#Loop over cycles to get structures:
#Applying Baker transformation:
mnum =[]
#open(fpdb,'w').close()
for i in (os.listdir(spath)):
    pdb =   PDBFile(os.path.join(spath,i))
    ENE = energy.energy(pdb)
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
            energy.outp(pdb,mname,mno,mpath,path)
            break
    

        print('----------------------------------------------------')
    
    if ENE<thresh:
        
        qmc.accept_move(pdb,mname,mno,mpath,path) 

       # print('****************************************************')
#Format conversions:
flibs.pdbread()
#Cluster conformations:
flibs.clust()
flibs.dihs()



#qgen.writeang(s,npar,b)





job_end = datetime.now()
print("Job Ended at:",job_end)
jend = process_time()

print("Time Elapsed (s):",(jend-jbegin))
