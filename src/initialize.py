#Intilalize input variables:

import os
import numpy as np





class VARINIT:
    def read_inp(self,inpfile):

        with open(inpfile,'r')as f:
            inp=[]
            for i,j in enumerate(f):
                if i==0:
                    self.path = str(j.split()[0])
                if i==1:
                    self.mname = str(j.split()[0])
                if i==3:
                    self.nstr=int(j.split()[0])
                if i==10:
                    self.indq=int(j.split()[1])
                if i==11:
                    self.npars=int(j.split()[1])
                if i==12:
                    self.search_sp=float(j.split()[1])
          
                

    def molsread(self,inpfile,nstr,spath):
        with open(inpfile) as f:
            for i,j in enumerate(f):
                lines=i
        sp = int(np.round(lines/nstr))
        #print(sp)
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


class STRPROP:

    def readf(self,pdb):
        with open(pdb,r) as ff:
            for i in ff.readlines():
                if i.startswith('MODEL'):
                    self.MNO = int(i.split()[1])
                if i.startwith('ATOM'):
                    self.atom_name.append(str(i.split()[2]))
                    self.atom_no.append(int(i.split()[1]))
                    self.resname.append(str(i.split()[3]))
                    self.resno.append(int(i.split()[4]))
                    self.x.append(float(i.split()[5]))
                    self.y.append(float(i.split()[6]))
                    self.z.appendz(float(i.split()[7]))
    
        
         
