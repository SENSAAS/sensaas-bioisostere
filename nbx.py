#!/usr/bin/python3.7

import os, sys
import re

ARGV=[]
try:
    ARGV.append(sys.argv[1])
    filesdf=ARGV[0]
except:
    print('usage: .py <sdf file to read>\n')
    quit()

# sys.argv[0] is the name of the program itself
filesdf=sys.argv[1]

#######################################
#def (filesdf) count number of X atoms
def nbx(fsdf):
    #read input file
    tabsdf=open(fsdf,'r')
    getstr=tabsdf.read().split('\n')
    tabsdf.close()
    tabLignesSdf=[]
    compt=3
    while(compt < len(getstr)):
        tabLignesSdf.append(re.split('\s+', getstr[compt].strip()))
        compt=compt+1
    testspace=[]
    testspace.append(re.split('', getstr[3]))
    if(len(tabLignesSdf[0][0]) > 2):#attached
        if(testspace[0][1]==' '):
            tabLignesSdf[0][1]=tabLignesSdf[0][0][2:]
            tabLignesSdf[0][0]=tabLignesSdf[0][0][0:2]
        elif(testspace[0][4]!=' '):
            tabLignesSdf[0][1]=tabLignesSdf[0][0][3:]
            tabLignesSdf[0][0]=tabLignesSdf[0][0][:3]
        nbatom=int(tabLignesSdf[0][0])
        nbbond=int(tabLignesSdf[0][1])
    else:
        nbatom=int(tabLignesSdf[0][0])
        nbbond=int(tabLignesSdf[0][1])
    #print("nbatom= %3s nbbond= %3s" % (nbatom,nbbond))
    tabR= {'Cl':1.75, 'Br':1.85, 'I':1.98, 'F':1.47, 'H':'%.2f' % 1.20, 'X':'%.2f' % 1.10}
    #Extract coordinates, atom type
    getA=[]
    getA.append('')
    nbxatoms=0
    compt=1
    #arr_xyz=np.empty(shape=[nbatom,3], dtype='float64')
    while (compt <= nbatom):
        #arr_xyz[compt-1,0]=float(tabLignesSdf[compt][0])
        #arr_xyz[compt-1,1]=float(tabLignesSdf[compt][1])
        #arr_xyz[compt-1,2]=float(tabLignesSdf[compt][2])
        getA.append(tabLignesSdf[compt][3])
        if(getA[compt] in tabR):
            if(getA[compt]=='X'):
                nbxatoms=nbxatoms+1
        compt=compt+1
    return nbxatoms

############################################
#MAIN
nx=nbx(filesdf)
print("%s" % nx)


