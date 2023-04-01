import os
from initialize import VARINIT

vinit = VARINIT()




vinit.read_inp('usermols.inp')



#Set files and folders:
class PATHS:
    def fpaths(self):
        mode = 0o777
        sampled = 'QMC_SAMPLED'
        minimized = 'QMC_MINIMIZED'
        self.spath = os.path.join(vinit.path,sampled)
        self.mpath = os.path.join(vinit.path,minimized)
        os.makedirs(self.spath,mode,exist_ok=True)
        os.makedirs(self.mpath,mode,exist_ok=True)
        #For accepted conformations:
        ffpdb = "{}_sobol.out".format(vinit.mname)
        fpdb = os.path.join(vinit.path,ffpdb) 
        self.molsfile = '{}_mols.pdb'.format(vinit.mname)

