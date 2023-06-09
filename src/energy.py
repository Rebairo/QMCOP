# Imports PDB conformation file and calculates energy:
import pandas as pd
from openmm.app import *
from openmm import *
from openmm.unit import *
import qgen
import numpy as np
import os
#from pdbfixer import PDBFixer




class ENERGY:

    def molsread(self,inpfile):
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

        
    def qmcene(self,pdb):
    
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


    def outp(self,pdb,mname,mno,mpath,path):

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

