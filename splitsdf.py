#!/usr/bin/python3.7

import os, sys
import re


ARGV=[]
try:
    ARGV.append(sys.argv[1])
    filesdf=ARGV[0]
except:
    print('usage: .py <sdf file to separate> <generic name>\n')
    quit()

# sys.argv[0] is the name of the program itself
filesdf=sys.argv[1]
genericname=sys.argv[2]

#######################################
#def (filesdf,generic-outputname) extract sdf no i into outputname i
def splitsdf(fsdf,osdf):
    #read input file
    tabsdf=open(fsdf,'r')
    getstr=tabsdf.read().split('\n')
    tabsdf.close()
    tabLignesSdf=[]
    whichend="$$$$"
    nbmol=0
    compt=0
    while(compt < len(getstr)):
        ri=nbmol+1
        namei=osdf + "_" + str(ri) + ".sdf"
        tabLignesSdf.append(getstr[compt])
        if(whichend in getstr[compt]):
            nbmol=nbmol+1
            ofile=open(namei, 'w')
            li=0
            while(li < len(tabLignesSdf)):
                ofile.write("%s\n" % tabLignesSdf[li])
                li=li+1
            ofile.close()
            tabLignesSdf=[]
        compt=compt+1
    if(nbmol==0):
        #print .mol file
        ri=1
        namei=osdf + "_" + str(ri) + ".sdf"
        li=0
        while(li < len(tabLignesSdf)):
            ofile.write("%s\n" % tabLignesSdf[li])
            li=li+1
    ofile.close()
    return

############################################
#MAIN
splitsdf(filesdf,genericname)

