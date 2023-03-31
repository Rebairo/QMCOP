# Imports PDB conformation file and calculates energy:
import pandas as pd
from openmm.app import *
from openmm import *
from openmm.unit import *
import qgen
import numpy as np
import os
#from pdbfixer import PDBFixer


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

        
def energy(pdb):
    
    forcefield = ForceField('amber96.xml')
#    forcefield = ForceField('amber14-all.xml','amber10_obc.xml')
    
#    forcefield = ForceField('amber14-all.xml','tip3p.xml')
    modeller = Modeller(pdb.topology,pdb.positions)
    modeller.addHydrogens(forcefield)
#    for ii in modeller.positions:
#            ii = (ii).in_units_of(nanometers)
#    print(modeller.positions) 

#        modeller.addSolvent(forcefield,padding=1.0*nanometers)
#    modeller.addSolvent(forcefield,numAdded=50)
#    modeller.addSolvent(forcefield,model='tip3p')
        #modeller.addMissingAtoms(forcefield)
    system = forcefield.createSystem(modeller.topology,ignoreExternalBonds=True)
    integrator = LangevinMiddleIntegrator(300*kelvin, 5/picosecond, 0.004*picoseconds)
    simulation = Simulation(modeller.topology, system,integrator)
    simulation.context.setPositions(modeller.positions)
    state = simulation.context.getState(getEnergy=True)
        #o = int(i/sp)
    ENE = state.getPotentialEnergy()
    ene = str(ENE)
    ene =  ene.split()[0]
    ene = ene[:8]
    ene = (float(ene)/10)
    ene = float("{:.3f}".format(ene))
#        with open(score,'a')as f:
#            f.write('%s %s %s %s\n'%("Model:",i, "Energy:",state.getPotentialEnergy()))
       
#        print('%s %s %s'%('QMC_SAMPLED_Structure:',i,state.getPotentialEnergy()))
    return ene


def outp(pdb,mname,mno,mpath,path):

#    forcefield = ForceField('amber10.xml','amber10_obc.xml')
    forcefield = ForceField('amber96.xml')
#    forcefield = ForceField('amber14-all.xml','tip3p.xml')
    modeller = Modeller(pdb.topology,pdb.positions)
    modeller.addHydrogens(forcefield)
#    modeller.addSolvent(forcefield,numAdded=50)
#    modeller.addSolvent(forcefield,model='tip3p')
    system = forcefield.createSystem(modeller.topology,ignoreExternalBonds=True)
    integrator = LangevinMiddleIntegrator(300*kelvin, 5/picosecond, 0.004*picoseconds)
    simulation = Simulation(modeller.topology, system,integrator)
    simulation.context.setPositions(modeller.positions)
    simulation.minimizeEnergy()
    positions = simulation.context.getState(getPositions=True).getPositions()
    output = os.path.join(mpath,'mini_{}.pdb'.format(mno))
    PDBFile.writeFile(simulation.topology, positions, open(output, 'w'))
    state = simulation.context.getState(getEnergy=True)
    ENE = state.getPotentialEnergy()
    ene = str(ENE)
    ene =  ene.split()[0]
    ene = ene[:8]
    ene = (float(ene)/10)
    ene = float("{:.3f}".format(ene))
    mout = "minimiz.out"
    minout = os.path.join(path,mout)
    mpdb = "{}_tmp_mini.pdb".format(mname)
    minpdb = os.path.join(path,mpdb)

    with open(minpdb,'a') as q:
        q.write('%s %d %f\n'%('MODEL',mno,ene))
    PDBFile.writeFile(simulation.topology,positions,open(minpdb,'a'))
    with open(minout,'a') as f:
        f.write('%s %s %s\n'%('QMC_MINIMIZED_Structure:',mno,ene))


# def minimiz(spath,mname):
    # for i in os.listdir(spath):

       # # pdb = PDBFile(os.path.join(path,i))
        # pdb = PDBFile(os.path.join(spath,i))
# #        forcefield = ForceField('amber14-all.xml','amber14/tip3p.xml')
        # forcefield = ForceField('amber14-all.xml')
        # modeller = Modeller(pdb.topology,pdb.positions)
        # modeller.addHydrogens(forcefield)
# #        modeller.addSolvent(forcefield,padding=1.0*nanometers)
# #        modeller.addSolvent(forcefield,numAdded=5)
# #        modeller.addSolvent(forcefield,model='tip3p')
       # #modeller.addMissingAtoms(forcefield)
        # system = forcefield.createSystem(modeller.topology,ignoreExternalBonds=True)
        # integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.004*picoseconds)
        # simulation = Simulation(modeller.topology, system,integrator)
        # simulation.context.setPositions(modeller.positions)
        # simulation.minimizeEnergy()
        # positions = simulation.context.getState(getPositions=True).getPositions()
        # #o = (i/sp)
        # output = os.path.join(mpath,'mini_{}'.format(i))
        # PDBFile.writeFile(simulation.topology, positions, open(output, 'w'))
        # state = simulation.context.getState(getEnergy=True)
        # ENE = state.getPotentialEnergy()
        # ene = str(ENE)
        # ene =  ene.split()[0]
        # ene = ene[:8]
        # ene = (float(ene)/10)
        # ene = float("{:.3f}".format(ene))
        # mout = "minimiz.out"
        # minout = os.path.join(xpath,mout)
        # mpdb = "{}_tmp_mini.pdb".format(mname)
        # minpdb = os.path.join(xpath,mpdb)
        # mnumber = i.partition('.')[0] 
        # mnumber = int(mnumber)

        # with open(minpdb,'a') as q:
            # q.write('%s %d %f\n'%('MODEL',mnumber,ene))
        # PDBFile.writeFile(simulation.topology,positions,open(minpdb,'a'))
        # with open(minout,'a') as f:
   

        #print('%s %s %s'%('QMC_MINIMIZED_Structure:',i,ENE))






#os.makedirs(spath,mode,exist_ok=True)
#os.makedirs(mpath,mode,exist_ok=True)
#os.mkdir(spath,mode)
#os.mkdir(mpath,mode)

#pdb =   PDBFile('check.pdb')
#K = energy(pdb)
#print(K)
#sp = molsread(os.path.join(xpath,molsfile))
#fix_pdb(path)
#print('\n')
#print('-------------------------BEGIN SAMPLING ---------------------\n') 
#print('\n')
#qmcscore(spath)
#print('\n')
#print('-------------------------SAMPLING DONE---------------------\n') 
# print('-------------------------BEGIN MINIMIZATION---------------------\n') 
# print('\n')
# minimiz(spath,mname)
# print('\n')
# print('-------------------------MINIMIZATION DONE---------------------\n') 
# #print(sp)
