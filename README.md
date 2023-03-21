# QMCOP
The Repo contains code that predicts polypeptide conformations from their sequence:

1. src/ - contains the source code.  
2. res/ - contains the results for corresponding polypeptides.    

In results folder:  
    # Energy is in (1/10) of (KJ/mol). 
    1. *_mini.pdb - contains final minimized conformations.  
    2. *_clut.pdb - contains clusters and their nearest structures.  
    3. *_pca.dat - contains PCA components of Dihedral angles [Energy,PCA1,PCA2].  
    4. *_eed.dat - contains Energy, end-end dist, Radgyr, Sq.Radgyr.  

# Instructions:  
1. Requirements: OpenMM library, Intel Fortran Compiler.  
2. Compile libraries ($sh make2),($sh make3),($sh make4),($sh make5)  
3. Input file - "src/usermols.inp"

#Input file description:  
File : usermols.inp.  
Lines:  
1 - location of the results directory.  
2 - Molecule name.  
3 - Peptide sequence.  
4 - No of conformations to be generated. (should be in 2^k). 
11 - ntors (No of torsions to sample).  
12 - search space (default 0-360). 

# Execution:  
To execute the program ($python main.py). 
