import os
import sys




# Fortran libs to be called:

class FORTLIBS:
    def conf(self):
        os.system("./lmols")

    def pdbread(self):    
        os.system("./pdbread")

    def clust(self):
        os.system("./clust")

    def dihs(self):
        os.system("./diheds")
