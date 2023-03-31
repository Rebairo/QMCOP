#Intilalize input variables:

import os
import numpy as np


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

def molsread(inpfile,nstr,spath):
    with open(inpfile) as f:
        for i,j in enumerate(f):
            lines=i
    sp = int(np.round(lines/nstr))
    print(sp)
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

