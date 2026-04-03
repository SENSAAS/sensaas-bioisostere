#!/usr/bin/python3.7

import os, shutil
import os, sys
import os.path
from sys import platform
import re

ARGV=[]
try:
    ARGV.append(sys.argv[1])
    filesdf=ARGV[0]
except:
    print('usage: .py <ordered-filtered-catsensaas.sdf> <name of fragments to combine eg: fgt_ (or part_ for part_*.sdf) > <file name of the part to replace (eg: part_3.sdf or fgt_1.sdf)')
    quit()

# sys.argv[0] is the name of the program itself
sensaasexe=sys.argv[0]
sensaasexe=re.sub('\/build-bioisostere\.py','',sensaasexe)

whichexe='linux'
if(whichexe in platform):
    # linux
    leaexec=sensaasexe + "/lea3d-combination.pl"
elif platform == "darwin":
    leaexec=sensaasexe + "/lea3d-combination.pl"
else:
    #windows with conda
    leaexec="perl " + sensaasexe + "\lea3d-combination.pl"

#print("perl %s " % (leaexec))

catsensaas=sys.argv[1]
genericname=sys.argv[2]
fgt=sys.argv[3]

#######################################
def isexist(file_path):
    return os.path.exists(file_path)

def isexistnotempty(file_path):
    if(os.path.exists(file_path)):
        if(os.path.getsize(file_path)==0):
            fempty=0
        else:
            fempty=1
    else:
        fempty=0
    return fempty

#######################################
#def (filesdf) output=number of molecules
def nbsdf(fsdf):
    tabsdf=open(fsdf,'r')
    getstr=tabsdf.read().split('\n')
    tabsdf.close()
    whichend="$$$$"
    nbmol=0
    compt=0
    while(compt < len(getstr)):
        if(whichend in getstr[compt]):
            nbmol=nbmol+1
        compt=compt+1
    if(nbmol==0): #if .mol format
        nbmol=1
    return nbmol

#######################################
#def (filesdf,i,outputname) extract sdf no i from a sdf file into outputname
def searchsdfi(fsdf,fi,osdf):
    #initiate output file
    ofile=open(osdf, 'w')
    #read input file
    tabsdf=open(fsdf,'r')
    getstr=tabsdf.read().split('\n')
    tabsdf.close()
    tabLignesSdf=[]
    whichend="$$$$"
    nbmol=0
    compt=0
    while(compt < len(getstr)):
        tabLignesSdf.append(getstr[compt])
        if(whichend in getstr[compt]):
            nbmol=nbmol+1
            if(nbmol==fi):
                #print
                li=0
                while(li < len(tabLignesSdf)):
                    ofile.write("%s\n" % tabLignesSdf[li])
                    li=li+1
            tabLignesSdf=[]
        compt=compt+1
    if(nbmol==0):
        #print .mol file
        li=0
        while(li < len(tabLignesSdf)):
            ofile.write("%s\n" % tabLignesSdf[li])
            li=li+1
    ofile.close()
    return

#######################################
# MAIN program

nbinput=nbsdf(catsensaas)
if(isexist("build-bioisostere.sdf")):
        os.remove("build-bioisostere.sdf")
shutil.copyfile(fgt,"copyfgt.sdf")
os.remove(fgt)
output="build-bioisostere.sdf"
mfile=open(output, 'w')
i=1
while(i <= nbinput):
    searchsdfi(catsensaas,i,fgt)
    if(isexist("combination.sdf")):
        os.remove("combination.sdf")
    cmd = '%s %s %s 2 0 ' % (leaexec,fgt,genericname)
    os.system(cmd)
    if(isexistnotempty("combination.sdf")):
        solfile=open('combination.sdf', 'r')
        stran=solfile.readlines()
        solfile.close()
        for f in stran:
            mfile.write(f)
    i=i+1

mfile.close()
nboutput=nbsdf("build-bioisostere.sdf")
print("Output in build-bioisostere.sdf %s molecules" % nboutput)

