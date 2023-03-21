import os
import qgen





nstr,indq,npars,xpath,mname = qgen.read_inp('usermols.inp')
minifile = '{}_mini.pdb'.format(mname)
cd = os.getcwd()
mpdb = os.path.join(xpath,minifile)
zz=[]
with open(mpdb) as f:
    for i,j in enumerate(f):
        if j.split()[0]=='REMARK':
            zz.append(i)
    lines = f.readlines()
print(zz)

with open(mpdb) as h:

    for ix in zz:
        print(h.seek(ix))


    
