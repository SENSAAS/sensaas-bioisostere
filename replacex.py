#!/usr/bin/python3.7

import os, sys
import re

ARGV=[]
try:
    ARGV.append(sys.argv[1])
    filesdf=ARGV[0]
except:
    print('usage: .py <sdf file to read and replace X by H>\n')
    quit()

# sys.argv[0] is the name of the program itself
filesdf=sys.argv[1]

#######################################
#def (filesdf) count number of X atoms
def replacex(fsdf):
    #read input file
    tabsdf=open(fsdf,'r')
    getstr=tabsdf.read().split('\n')
    tabsdf.close()
    tabR= {'Cl':1.75, 'Br':1.85, 'I':1.98, 'F':1.47, 'H':'%.2f' % 1.20, 'X':'%.2f' % 1.10}
    ofile=open("replacex.sdf", 'w')
    whichend="$$$$"
    li=0
    nbmol=0
    new=1
    tabLignesSdfbis=[]
    while(li < len(getstr)):
        if(new==1 and (li+3) < len(getstr)):
            point=li
            #print("li %s" % li)
            #print("li %s" % getstr[li])
            new=0
            nbmol=nbmol+1
            tabLignesSdf=[]
            ligneab=li+3
            tabLignesSdf.append(re.split('\s+', getstr[ligneab].strip()))
            testspace=[]
            testspace.append(re.split('', getstr[ligneab]))
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
        #check atom type
        tabLignesSdfbis.append(re.split('\s+', getstr[li].strip()))
        if(new==0 and li  > (point+3) and li <= (point+nbatom+3)):
            #print ("point %s %s" % (point,getstr[li]))
            elt=tabLignesSdfbis[li][3]
            #print ("elt %s" % elt)
            if(elt=="X"):
                #print ("ici %s" % getstr[li])
                #replace by H
                line=""
                tabLine=re.split(' +', getstr[li].strip())
                line=line + str('%10s%10s%10s' % (tabLine[0],tabLine[1],tabLine[2]))
                tabLine[3]="H"
                line=line + ' ' + str('%1s' % tabLine[3]) + ' '
                maxi=len(tabLine)
                for i in range(4,maxi):
                    line=line + ' ' + str('%2s' % tabLine[i])
                ofile.write(line+'\n')
                #print("ici %s\n" % line)
            else:
                ofile.write("%s\n" % getstr[li])
        elif(new==0):
            if(whichend in getstr[li]):
                ofile.write("$$$$\n")
                new=1
            else:
                ofile.write("%s\n" % getstr[li])
        li=li+1
    ofile.close()
    #print("nb mol %s" % nbmol)
    return 

############################################
#MAIN
replacex(filesdf)

