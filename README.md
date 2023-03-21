# QuPoly
The Repo contains code that predicts polypeptide conformations from their sequence:

1. src/ - contains the source code.  
2. res/ - contains the results for corresponding polypeptides.    

In results folder:  
    # Energy is in (1/10) of (KJ/mol)
    1. *_mini.pdb - contains final minimized conformations.  
    2. *_clut.pdb - contains clusters and their nearest structures.  
    3. *_pca.dat - contains PCA components of Dihedral angles [Energy,PCA1,PCA2].  
    4. *_eed.dat - contains Energy, end-end dist, Radgyr, Sq.Radgyr.  
