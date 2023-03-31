import os
import initialize










#Set files and folders:

def path(path,mname):
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
    molsfile = '{}_mols.pdb'.format(mname)

    return spath,mpath,molsfile
