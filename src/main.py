# MAIN QMC PROGRAM:
import sys
import qgen
import scipy
import os
import numpy as np
from scipy.stats import qmc
import math
import pandas as pd
import random
import energy

from openmm.app import *
from openmm import *
from openmm.unit import *
#Provide an index to choose Sobol-1/RandomLHS-2/OrhtogonalLHS-3/Halton-4 quasi-random sequences:








#Function defs:
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

def molsread(inpfile):
    with open(inpfile) as f:
        for i,j in enumerate(f):
            lines=i
    sp = int(np.round(lines/nstr))

    spdb=None
    with open(inpfile) as f:

        for i,j in enumerate(f):
            if i%sp==0:
                if spdb:
                    spdb.close()
                o = int(i/sp)
                s_pdb=os.path.join(spath,'{}.pdb'.format(o))
                spdb = open(s_pdb,'w')
            spdb.write(j)
        if spdb:
            spdb.close()
    return sp


#Get Input Vars:
nstr,ind,npar,path,mname,search_sp  =  read_inp('usermols.inp')

#--------------------------------#
#Threshold:
thresh = 200.0
thresh2 = 1000.0


#Run Initial Sample:
if ind==1:
    s,b=qgen.sobol(npar,nstr)

#Initial Conformation generation:
os.system("./lmols")

#Set files and folders:
mode = 0o777
sampled = 'QMC_SAMPLED'
minimized = 'QMC_MINIMIZED'
spath = os.path.join(path,sampled)
mpath = os.path.join(path,minimized)
os.makedirs(spath,mode,exist_ok=True)
os.makedirs(mpath,mode,exist_ok=True)
#For accepted conformations:
ffpdb = "{}_sobol.out".format(mname)

fpdb = os.path.join(path,ffpdb)
#Use conformations and split them:
molsfile = '{}_mols.pdb'.format(mname)
sp = molsread(os.path.join(path,molsfile))
#Loop over cycles to get structures:
#Applying Baker transformation:
mnum =[]
open(fpdb,'w').close()
for i in (os.listdir(spath)):
    pdb =   PDBFile(os.path.join(spath,i))
    ENE = energy.energy(pdb)
    mno = int(i.partition('.')[0])



    print(ENE)
    bb=0
    if ENE>thresh2:
        while ENE>thresh2:
            bb = bb+1
        
            print('RESAMPLING MOVES \n')
            print('MODEL NO:',mno)
            print('CYCLE:',bb)
            mnum.append(mno)
        
            with open ('mols_q.out','r+') as ff:
                for j in ff.readlines():
                    if(mno==int(j.split()[0])):
                            hh = j
        #print(hh) 
            jj = hh.split()[1:]
        #udat = np.around(np.asarray(jj,dtype=float),decimals=2)
        #U = (np.random.random(npar))*(search_sp)
            U2,b = qgen.sobol(npar,1) 
        #print(U2)
        #sdat =(((udat-(1-2*U))%2))
            s2dat = np.transpose(U2)
            s2dat2 = s2dat.reshape((1,npar))
        #print(sdat)
        #print(sdat2)
            df2_udat = pd.DataFrame(s2dat2)
        #print(df_udat)
        #exit()
            df2_udat.to_csv("mols_q.out",sep="\t",header=None,index=True)
            os.system("./lmols")
        
        #sp = molsread(os.path.join(path,molsfile))
            pdb =   PDBFile(os.path.join(path,molsfile))
            ENE = energy.energy(pdb)
            print('ENERGY AFTER RESAMPLING:',ENE)
            print('--------------------------')


    cc=0
    if ENE>thresh and ENE<thresh2:
        while ENE>thresh:
            cc = cc +1
            print('APPLYING MOD M SHIFT TRANSFORMATION')
            mnum.append(mno)
        
            with open ('mols_q.out','r+') as ff:
                for j in ff.readlines():
                    if(mno==int(j.split()[0])):
                            hh = j
        #print(hh) 
            jj = hh.split()[1:]
            udat = np.around(np.asarray(jj,dtype=float),decimals=2)
            U = (np.random.random(npar))*(search_sp)
        #print(U)
            sdat =(((udat-(U))))
        #sdat = np.transpose(sdat)
            sdat2 = sdat.reshape((1,npar))
        #print(sdat)
        #print(sdat2)
            df_udat = pd.DataFrame(sdat2)
        #print(df_udat)
        #exit()
            df_udat.to_csv("mols_q.out",sep="\t",header=None,index=True)
            os.system("./lmols")
        
        #sp = molsread(os.path.join(path,molsfile))
            pdb =   PDBFile(os.path.join(path,molsfile))
            ENE = energy.energy(pdb)
            print('ENERGY AFTER MOD SHIFT TRANSFORMATION:',ENE)
            print("Cycle:",cc)
            
    if ENE<thresh:
        print('MINIMISING OPTIMAL CONFORMATION')
        energy.outp(pdb,mname,mno,mpath,path)
        print('OPTIMAL CONFORMATION GENERATED')
            






# Shift unfavourable states by Bakers map:
#sdat = np.zeros([len(mnum),npar])
# SHift:
#U = ((np.round(np.random.rand(npar),decimals=2))*(search_sp/3))
#sdat = np.transpose(udat)-np.vstack(U)
#sdat = np.transpose(sdat)

#print((sdat-udat))
#Append MODEL NO with SHIFTED DATA:
#mnum = (np.asarray(mnum,dtype=int))
#usdat = np.around(np.column_stack([(mnum),sdat]),decimals=2)

#df_usdat = pd.DataFrame(usdat)
#df_usdat.to_csv("mols_q.out",sep="\t",header=None,index=None)

#Recalculate ENERGY:

#print(np.transpose(sdat))
     
os.system("./pdbread")
os.system("./clust")
os.system("./diheds")
#qgen.writeang(s,npar,b)
