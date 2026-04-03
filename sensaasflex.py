#!/usr/bin/python3.7

import os, shutil
import os, sys
import os.path
import re
import numpy as np
from sys import platform
import random
import math
from SDForPDBtoDots import *
import TemplateS

ARGV=[]
try:
    ARGV.append(sys.argv[1])
    filesdf=ARGV[0]
except:
    print('usage: sensaasflex.py target-sdffile-name source-sdffile-name \neg:\nsensaasflex.py P04035-7.sdf P04035-7-confs1.sdf ')
    quit()

# sys.argv[0] is the name of the program itself
sensaasexe=sys.argv[0]
sensaasexe=re.sub('sensaasflex\.py','',sensaasexe)

whichexe='linux'
if(whichexe in platform):
    # linux
    rdkitminexe=sensaasexe + "rdkit-min-iter.py"
    sensaasexecicp=sensaasexe + "sensaasCICP.py"
    sensaasexegcicp=sensaasexe + "sensaas.py"
    leaexec=sensaasexe + "lea3d-combination.pl"
    leaexefgtagg=sensaasexe + "lea3d-MAKE_FGTS_AGGREG.pl"
elif platform == "darwin":
    # OS X - linux version?
    rdkitminexe=sensaasexe + "rdkit-min-iter.py"
    sensaasexecicp=sensaasexe + "sensaasCICP.py"
    sensaasexegcicp=sensaasexe + "sensaas.py"
    leaexec=sensaasexe + "lea3d-combination.pl"
    leaexefgtagg=sensaasexe + "lea3d-MAKE_FGTS_AGGREG.pl"
else:
    #windows with conda
    rdkitminexe="python " + sensaasexe + "rdkit-min-iter.py"
    sensaasexecicp="python " + sensaasexe + "sensaasCICP.py"
    sensaasexegcicp="python " + sensaasexe + "sensaas.py"
    leaexec="perl " + sensaasexe + "lea3d-combination.pl"
    leaexefgtagg="perl " + sensaasexe + "lea3d-MAKE_FGTS_AGGREG.pl"

#print("exe %s and %s and %s" % (sensaasexe,sensaasexegcicp,sensaasexecicp))
#print("perl %s and %s " % (leaexec,leaexefgtagg))
#print("exe %s" % (rdkitminexe))

target=sys.argv[1]
source=sys.argv[2]

nbarg=len(sys.argv)
if(nbarg==4):
    targetdots=sys.argv[3]
    threshold=0.3
elif(nbarg==5):
    targetdots=sys.argv[3]
    threshold=float(sys.argv[4])
else:
    targetdots=""
    threshold=0.3

nbDots=294
#nbDots=182
TemplateDots=np.empty(shape=[nbDots,3], dtype='float64')
tabo=[]
tabo=TemplateS.tab
TemplateDots=np.array(tabo)

#number of rounds (usually = 2)
nbrounds=2
# static=1 (not reorientated) or static=0 (to rotate and translate source structure)
static=0
flagcleanfiles=1
verbose=0

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
#def randomizesdf(inputfile,outputfile)
def randomizesdf(fsdf,osdf):
    #create random matrix 4D with 16 values
    alpha=random.randint(0,180)
    a1=math.cos(alpha)
    a2=math.sin(alpha)
    a3=0
    b1=-1*(math.sin(alpha))
    b2=math.cos(alpha)
    b3=0
    c1=0
    c2=0
    c3=1
    a4=random.uniform(0,5)
    t2=random.uniform(0,1)
    if (t2<=0.5):
        a4=a4*-1
    b4=random.uniform(0,5)
    t2=random.uniform(0,1)
    if (t2<=0.5):
        b4=b4*-1
    c4=random.uniform(0,5)
    t2=random.uniform(0,1)
    if (t2<=0.5):
        c4=c4*-1
    d1=0
    d2=0
    d3=0
    d4=1
    #read input file
    tabsdf=open(fsdf,'r')
    getstr=tabsdf.read().split('\n')
    tabsdf.close()
    tabLignesSdf=[]
    tabLignesSdf.append(re.split('\s+', getstr[3].strip()))
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
    #initiate output file
    ofile=open(osdf, 'w')
    compt=0
    while(compt < len(getstr)):
        if(compt > 3 and compt <= nbatom+3):
            tabLine=re.split(' +', getstr[compt].strip())
            x=float(tabLine[0])
            y=float(tabLine[1])
            z=float(tabLine[2])
            x1 = x*a1 + y*a2 + z*a3 + a4
            y1 = x*b1 + y*b2 + z*b3 + b4
            z1 = x*c1 + y*c2 + z*c3 + c4
            line=""
            tabLine[0]=str('%5.3f' % (x1))
            tabLine[1]=str('%5.3f' % (y1))
            tabLine[2]=str('%5.3f' % (z1))
            line=line + str('%10s%10s%10s' % (tabLine[0],tabLine[1],tabLine[2]))
            elt=tabLine[3]
            longa=len(tabLine[3])
            if(longa==1):
                line=line + ' ' + str('%1s' % tabLine[3]) + ' '
            else:
                line=line + ' ' + str('%2s' % tabLine[3])
            maxi=len(tabLine)
            for i in range(4,maxi):
                line=line + ' ' + str('%2s' % tabLine[i])
            ofile.write(line+'\n')
        else:
            if(compt==(len(getstr)-1)):
                ofile.write(getstr[compt])
            else:
                ofile.write(getstr[compt]+'\n')
        compt=compt+1
    ofile.close()
    return

#######################################
#def xyzsdf(inputfile)
def xyzsdf(fsdf):
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
    nbsmall=0
    nblarge=0
    compt=1
    arr_xyz=np.empty(shape=[nbatom,3], dtype='float64')
    while (compt <= nbatom):
        arr_xyz[compt-1,0]=float(tabLignesSdf[compt][0])
        arr_xyz[compt-1,1]=float(tabLignesSdf[compt][1])
        arr_xyz[compt-1,2]=float(tabLignesSdf[compt][2])
        getA.append(tabLignesSdf[compt][3])
        if(getA[compt] in tabR):
            nbsmall=nbsmall+1
        else:
            nblarge=nblarge+1
        compt=compt+1
    return arr_xyz, nbatom, nblarge

#######################################
#def read slog from sensaas output
def readscores(filelog):
    #read scores
    logfile=open(filelog, 'r')
    lignelog=logfile.readlines()
    logfile.close()
    score=lignelog[-1]
    tabscore=[]
    tabscore.append(re.split('\s+', score.strip()))
    gfh=float(tabscore[0][7])
    g=float(tabscore[0][1])
    h=float(tabscore[0][5])
    c=float(tabscore[0][3])
    return gfh,g,h,c

#######################################
#def scorerecombination(tmpsa file)
def scorerecombination(tmpsafile):
    gfithfitco=0
    gfitco=0
    hfitco=0
    cfitco=0
    if(isexist("combination.sdf")):
        os.remove("combination.sdf")
    if(isexist("combination-min.sdf")):
        os.remove("combination-min.sdf")
    #print("%s %s make_fgts_ 2 " % (leaexec,tmpsafile))
    cmd = '%s %s make_fgts_ 2 ' % (leaexec,tmpsafile)
    os.system(cmd)
    if(isexistnotempty("combination.sdf")):
        #clean up by few iteration in minimization (rdkitmin)
        cmd = '%s combination.sdf combination-min.sdf' % (rdkitminexe)
        os.system(cmd)
        if(isexistnotempty("combination-min.sdf")):
            cmd = '%s dot tmpt-dots.pdb sdf combination-min.sdf slog optim %s' % (sensaasexecicp,threshold)
            os.system(cmd)
            gfithfitco,gfitco,hfitco,cfitco=readscores("slog")
            #os.rename("Source_tran.sdf","Source_tran-conf.sdf")
            shutil.copyfile("Source_tran.sdf","Source_tran-conf.sdf")
            os.remove("Source_tran.sdf")
    return gfithfitco,gfitco,hfitco,cfitco

#######################################
#def dotspdb2fgts(targetptcloud,sourcefile,distm,output-name)
def dotspdb2fgts(targetptcloud,sourcefile,distm,oname):
    if(isexist("make_fgts_aggreg.sdf")):
        os.remove("make_fgts_aggreg.sdf")
    if(isexist("make_fgts.sdf")):
        os.remove("make_fgts.sdf")
    cmd = '%s %s ' % (leaexefgtagg,sourcefile)
    os.system(cmd)
    if(isexistnotempty("make_fgts_aggreg.sdf")):
        #os.rename("make_fgts_aggreg.sdf","make_fgts.sdf")
        shutil.copyfile("make_fgts_aggreg.sdf","make_fgts.sdf")
        os.remove("make_fgts_aggreg.sdf")
    nbfgtagg=nbsdf("make_fgts.sdf")

    if(nbfgtagg > 1):
        #read full target point cloud
        filemol=open(targetptcloud,'r')
        getstr=[]
        getstr=filemol.read().split('\n')
        filemol.close()
        tabLignesPdb=[]
        tabLignesPdb.append('')
        compt=1
        while (compt < len(getstr)):
            tabLignesPdb.append(re.split('\s+', getstr[compt].strip()))
            compt=compt+1

        fgti=1
        while (fgti <= nbfgtagg):
            #get coordinates of fgti
            namesdfi="make_fgts_" + str(fgti) + ".sdf";
            searchsdfi("make_fgts.sdf",fgti,namesdfi)
            coord,nbatoms,nbheavy=xyzsdf(namesdfi)
            if(nbheavy > 1):
                #only if fragment is > CH3, OH, NH3, halogen...
                #create fgti dots file
                namedotsi=oname + "_" + str(fgti) + "-dots.pdb";
                #if(isexist(namedotsp)):
                #    os.remove(namedotsp)
                solfile=open(namedotsi, 'w')
                compt=1
                while (compt < len(tabLignesPdb)):
                    if (tabLignesPdb[compt][0] == 'HETATM' or tabLignesPdb[compt][0] == 'ATOM'):
                        xAtome=float(tabLignesPdb[compt][5])
                        yAtome=float(tabLignesPdb[compt][6])
                        zAtome=float(tabLignesPdb[compt][7])
                        ligned=getstr[compt]
                        #dist
                        ci=0
                        vu=0
                        while (vu==0 and ci < nbatoms):
                            getdot= math.sqrt((xAtome - coord[ci,0])**2 + (yAtome - coord[ci,1])**2 + (zAtome - coord[ci,2])**2)
                            if(getdot <= float(distm)):
                                vu=1
                            ci=ci+1
                        if(vu==1):
                            solfile.write("%s \n" % ligned)
                    compt=compt+1
                solfile.close

            fgti=fgti+1
    return

#######################################
#def redefinepointcloud 
def redefinepointcloud(p,pt,filept):
    #exclude points close to other fragments than p
    getx2=[]
    getx2.append('')
    gety2=[]
    gety2.append('')
    getz2=[]
    getz2.append('')
    nbpt2=0
    li=1
    while (li <= pt):
        if(li != p):
            named="Sourceconf_" + str(li) + "-dots.pdb"
            if(isexist(named)):
                getstr=[]
                filemol=open(named,'r')
                getstr=filemol.read().split('\n')
                filemol.close()
                tabLignesPdb=[]
                tabLignesPdb.append('')
                compt=1
                while (compt < len(getstr)):
                    tabLignesPdb.append(re.split('\s+', getstr[compt].strip()))
                    compt=compt+1
                compt=1
                while (compt < len(tabLignesPdb)):
                    if (tabLignesPdb[compt][0] == 'HETATM' or tabLignesPdb[compt][0] == 'ATOM'):
                        xAtome=float(tabLignesPdb[compt][5])
                        yAtome=float(tabLignesPdb[compt][6])
                        zAtome=float(tabLignesPdb[compt][7])
                        getx2.append(xAtome)
                        gety2.append(yAtome)
                        getz2.append(zAtome)
                        nbpt2=nbpt2+1
                    compt=compt+1
        li=li+1

    namedotsp="Sourceconf_" + str(p) + "-dots.pdb";
    #if(isexist(namedotsp)):
    #    os.remove(namedotsp)
    solfile=open(namedotsp, 'w')
    
    #read full target point cloud
    filemol=open(filept,'r')
    getstr=[]
    getstr=filemol.read().split('\n')
    filemol.close()
    tabLignesPdb=[]
    tabLignesPdb.append('')
    compt=1
    while (compt < len(getstr)):
        tabLignesPdb.append(re.split('\s+', getstr[compt].strip()))
        compt=compt+1

    compt=1
    nbpt=0
    while (compt < len(tabLignesPdb)):
        if (tabLignesPdb[compt][0] == 'HETATM' or tabLignesPdb[compt][0] == 'ATOM'):
            xAtome=float(tabLignesPdb[compt][5])
            yAtome=float(tabLignesPdb[compt][6])
            zAtome=float(tabLignesPdb[compt][7])
            nbpt=nbpt+1
            li=1
            vu=0
            while (li <= nbpt2):
                if(getx2[li]==xAtome and gety2[li]==yAtome and getz2[li]==zAtome):
                    vu=1
                li=li+1
            if(vu==0):
                #getstr indice == tabLignesPdb indice:
                ligned=getstr[compt]
                #print(ligned)
                solfile.write("%s \n" % ligned)
        compt=compt+1
    solfile.close
    return

#######################################
# MAIN program

#works with only one target and source file:
i=1
j=1
searchsdfi(target,i,"tmptf.sdf")
searchsdfi(source,j,"tmpsf.sdf")
tmps="tmpsf.sdf"
if(static==0):
    randomizesdf(tmps,tmps)

#target will be a dot file
if(targetdots != "" and isexistnotempty(targetdots)):
    shutil.copyfile(targetdots,"tmpt-dots.pdb")
else:
    molsurface("sdf","tmptf.sdf",TemplateDots,nbDots,1)
    os.remove("dots.xyzrgb")
    os.remove("dotslabel4.xyzrgb")
    os.remove("dotslabel3.xyzrgb")
    os.remove("dotslabel2.xyzrgb")
    os.remove("dotslabel1.xyzrgb")
    shutil.copyfile("dots.pdb","tmpt-dots.pdb")
    os.remove("dots.pdb")
    #os.rename("dots.pdb","tmpt-dots.pdb")

#initiate alignment: optimize Globally + locally
gfithfit=np.empty(shape=[nbrounds+2,4], dtype='float64')
stran=np.empty(shape=[nbrounds+2,1], dtype='str')
gfh=0
gf=0
hf=0
cf=0
cmd = '%s dot tmpt-dots.pdb sdf %s slog optim %s' % (sensaasexegcicp,tmps,threshold)
os.system(cmd)
gfithfit[0,0],gfithfit[0,1],gfithfit[0,2],gfithfit[0,3]=readscores("slog")
stran=[]
solfile=open('Source_tran.sdf', 'r')
stran=solfile.readlines()
solfile.close()
initstran=[]
initstran=stran
best=0
nbfgt=0

#initiate output
output="slogflex"
mfile=open(output, 'w')
#mfile.close()
if(verbose==1):
    mfile.write("%s - %s sensaas scores: %s \n" % (target,tmps,gfithfit[0,0]))
    print ("%s - %s sensaas scores: %s" % (target,tmps,gfithfit[0,0]))
else:
    mfile.write("Initial gfit= %s cfit= %s hfit= %s gfit+hfit= %s \n" % (gfithfit[best,1],gfithfit[best,3],gfithfit[best,2],gfithfit[best,0])) 

#parameters : gfithfit ref < 1.8 ; dots-sdf dist max=[2-4] ; make_fgt_aggreg id nbheavy <=2
if(gfithfit[0,0] < 1.8 and gfithfit[0,0] > 0 and nbrounds > 0):

    #tmpsa is the sdf file of the conformer to optimize 
    tmpsa="tmpsa.sdf"
    #os.rename("Source_tran.sdf",tmpsa)
    shutil.copyfile("Source_tran.sdf",tmpsa)
    os.remove("Source_tran.sdf")

    #split fgts and associate dots of the target point cloud (distance= 4 A)
    #3D graph  make_fgts_".$i.".sdf
    #Point cloud Source_".$i."-dots.pdb
    distmax=4
    distmaxb=2
    dotspdb2fgts("tmpt-dots.pdb",tmpsa,distmax,"Source")

    nbfgt=nbsdf("make_fgts.sdf")
    if(nbfgt >= 5 and nbrounds <= 2):
        if(verbose==1):
            print("Nb fgts >= 5 ; nbrounds set to 3 instead of %s" % nbrounds)
            mfile.write("Nb fgts >= 5 ; nbrounds set to 3 instead of %s \n" % nbrounds)
        nbrounds=3

    if(nbfgt <= 1):
        if(verbose==1):
            print("only one fgt: conformer search is skipped")
            mfile.write("only one fgt: conformer search is skipped \n")
    else:
        #work with fragments
        if(verbose==1):
            print("Optimize fragments [1 - %s] make_fgts_i.sdf on Source_i-dots.pdb (up to distance= %s A)" % (nbfgt,distmax))
            mfile.write("Optimize fragments [1 - %s] make_fgts_i.sdf on Source_i-dots.pdb (up to distance= %s A) \n" % (nbfgt,distmax))
        
        nbri=1
        while(nbri <= nbrounds):
            if(verbose==1):
                print("Round %s" % nbri)
                mfile.write("Round %s \n" % nbri)

            if(nbri > 1):
                #new tmpsa file : split fgts and associate dots of the target point cloud (distance= 4 A)
                dotspdb2fgts("tmpt-dots.pdb",tmpsa,distmax,"Source")
                #may change depending how acyclics are segmented
                nbfgt=nbsdf("make_fgts.sdf")

            fgtinfo=[]
            fgtinfo.append('')
            fgtscore=[]
            fgtscore.append('')
            #evaluate each fgt
            flags2=0
            ri=1
            while (ri <= nbfgt):
                #get scores of each fragment
                name="make_fgts_" + str(ri) + ".sdf"
                named="Source_" + str(ri) + "-dots.pdb"
                fgtinfo.append(int(0))
                fgtscore.append(float(0))
                if(isexistnotempty(named)):
                    #sensaas evaluation
                    cmd = '%s dot %s sdf %s slog eval %s' % (sensaasexegcicp,named,name,threshold)
                    os.system(cmd)
                    fgtscore[ri],gf,hf,cf=readscores("slog")
                    if(fgtscore[ri] > 0.3 and fgtscore[ri] < 1.2):
                        fgtinfo[ri]=1
                    elif(fgtscore[ri] <= 0.3):
                        fgtinfo[ri]=2
                        flags2=1
                elif(isexist(named)):
                    if(verbose==1):
                        print("fgt %s %s empty" % (ri,named))
                        mfile.write("fgt %s %s empty \n" % (ri,named))
                    fgtinfo[ri]=2
                    flags2=1
                else:
                    if(verbose==1):
                        print("fgt %s %s not found  = small atom with one atom" % (ri,named))
                        mfile.write("fgt %s %s not found  = small atom with one atom \n" % (ri,named))
                ri=ri+1

            #rework point cloud of fragments not aligned at all
            if(flags2==1):
                #target dots and source sdf are appart
                #split fgts and associate dots of the target point cloud (distance= 2 A)
                if(verbose==1):
                    print("redefine point cloud: use %s in dotspdb2fgts with %s A" % (tmpsa,distmaxb))
                    mfile.write("redefine point cloud: use %s in dotspdb2fgts with %s A \n" % (tmpsa,distmaxb))
                dotspdb2fgts("tmpt-dots.pdb",tmpsa,distmaxb,"Sourceconf")
                #allow to enlarge point cloud of fgt $fgtsinfo[$i]==2 after removing well aligned dots (Sourceconf-$i-dots.pdb) with dist=2A
                rj=1
                while (rj <= nbfgt):
                    if(fgtinfo[rj]==2):
                        #put all target dots not associated with any fgt at distmaxb A into the fgt point cloud Sourceconf-$i-dots.pdb
                        if(verbose==1):
                            print("redefine point cloud of %s (point cloud tmpt-dots.pdb minus other Sourceconf_j-dots.pdb (j!=i)" % rj)
                            mfile.write("redefine point cloud of %s (point cloud tmpt-dots.pdb minus other Sourceconf_j-dots.pdb (j!=i) \n" % rj)
                        redefinepointcloud(rj,nbfgt,"tmpt-dots.pdb")
                    rj=rj+1

            #optimize fgt ri in function of the configuration
            if(verbose==1):
                print("Refine fragments alignment:")
                mfile.write("Refine fragments alignment:")

            recombin=0
            ri=1
            while (ri <= nbfgt):
                if(fgtinfo[ri] > 0):
                    
                    name="make_fgts_" + str(ri) + ".sdf"
                    if(fgtinfo[ri]==2):
                        named="Sourceconf_" + str(ri) + "-dots.pdb"
                    else:
                        named="Source_" + str(ri) + "-dots.pdb"
                    if(verbose==1):
                        print("fgt %s ( %s and %s ) score initial= %s" % (ri,name,named,fgtscore[ri]))
                        mfile.write("fgt %s ( %s and %s ) score initial= %s \n" % (ri,name,named,fgtscore[ri]))
                    
                    #globally + locally on non-aligned point clouds
                    if(fgtinfo[ri]==2 and isexistnotempty(named)):
                        gfithfit4=0
                        cmd = '%s dot %s sdf %s slog optim %s' % (sensaasexegcicp,named,name,threshold)
                        os.system(cmd)
                        gfithfit4,gf,hf,cf=readscores("slog")

                        if(gfithfit4 > fgtscore[ri]):
                            #alignement of ri is improved
                            recombin=1
                            if(verbose==1):
                                print("fgt %s replaced - score after new alignt = %s" % (ri,gfithfit4))
                                mfile.write("fgt %s replaced - score after new alignt = %s \n" % (ri,gfithfit4))
                            #os.rename("Source_tran.sdf",name)
                            shutil.copyfile("Source_tran.sdf",name)
                            os.remove("Source_tran.sdf")
                            fgtscore[ri]=float(gfithfit4)
                            #success
                        else:
                            if(verbose==1):
                                print("fgt %s keep old - score after new alignt = %s" % (ri,gfithfit4))
                                mfile.write("fgt %s keep old - score after new alignt = %s \n" % (ri,gfithfit4))

                    #locally only if initial score of fgt is acceptable (fgtinfo[ri] == 1)
                    if(fgtinfo[ri]==1 and isexistnotempty(named)):
                        gfithfit2=0
                        cmd = '%s dot %s sdf %s slog optim %s' % (sensaasexecicp,named,name,threshold)
                        os.system(cmd)
                        gfithfit2,gf,hf,cf=readscores("slog")
                        if(gfithfit2 > fgtscore[ri]):
                            #alignement of ri is improved
                            recombin=1
                            if(verbose==1):
                                print("fgt %s replaced - score after refinement = %s" % (ri,gfithfit2))
                                mfile.write("fgt %s replaced - score after refinement = %s \n" % (ri,gfithfit2))
                            #os.rename("Source_tran.sdf",name)
                            shutil.copyfile("Source_tran.sdf",name)
                            os.remove("Source_tran.sdf")
                            fgtscore[ri]=float(gfithfit2)
                        else:
                            if(verbose==1):
                                print("fgt %s keep old - score after refinement = %s " % (ri,gfithfit2))
                                mfile.write("fgt %s keep old - score after refinement = %s \n" % (ri,gfithfit2))
                ri=ri+1

            #recombination into entire molecule
            gfithfit3=0
            if(recombin==1):
                if(verbose==1):
                    print("combination of fragments after individual alignment and refine alignement")
                    mfile.write("combination of fragments after individual alignment and refine alignment \n")
                gfithfit3,gf,hf,cf=scorerecombination(tmpsa)
                if(gfithfit3 > gfithfit[best,0]):
                    gfithfit[nbri,0]=gfithfit3
                    gfithfit[nbri,1]=gf
                    gfithfit[nbri,2]=hf
                    gfithfit[nbri,3]=cf
                    best=nbri
                    if(verbose==1):
                       print("round %s success new score %s (ref= %s)" % (nbri,gfithfit3,gfithfit[best,0]))
                       mfile.write("round %s success new score %s (ref= %s) \n" % (nbri,gfithfit3,gfithfit[best,0]))
                    else:
                        mfile.write("Round %s gfit= %s cfit= %s hfit= %s gfit+hfit= %s \n" % (nbri,gfithfit[best,1],gfithfit[best,3],gfithfit[best,2],gfithfit[best,0]))

                    if(isexist("Source_tran-conf.sdf")):
                        solfile=open('Source_tran-conf.sdf', 'r')
                        stran=solfile.readlines()
                        solfile.close()
                else:
                    if(verbose==1):
                        print("round %s fail new score= %s (ref= %s)" % (nbri,gfithfit3,gfithfit[best,0]))
                        mfile.write("round %s fail new score= %s (ref= %s) \n" % (nbri,gfithfit3,gfithfit[best,0]))
                    else:
                        mfile.write("Round %s gfit= %s cfit= %s hfit= %s gfit+hfit= %s \n" % (nbri,gf,cf,hf,gfithfit3))
                     #best not changed

            #stop or continue
            if(gfithfit[best,0] > 1 and best > 0):
                #stop because score was enough improved with the conformer
                nbri=nbrounds
                if(verbose==1):
                    print("stop because well improved")
                    mfile.write("stop because well improved \n")
            else:
                #continue 
                if(gfithfit[nbri,0] > 0 and isexist("Source_tran-conf.sdf")):
                    #os.rename("Source_tran-conf.sdf",tmpsa)
                    shutil.copyfile("Source_tran-conf.sdf",tmpsa)
                    os.remove("Source_tran-conf.sdf")
                    if(verbose==1):
                        print("next round will try to improve current Source_tran-conf.sdf")
                        mfile.write("next round will try to improve current Source_tran-conf.sdf \n")
                else:
                    #improve best 
                    if(best > 0):
                        catfile=open(tmpsa,'w')
                        for f in stran:
                            catfile.write(f)
                        catfile.close()
                        if(verbose==1):
                            print("next round will try to improve best solution %s with score= %s" % (gfithfit[best,0])) 
                            mfile.write("next round will try to improve best solution %s with score= %s \n" % (gfithfit[best,0]))
                    else:
                        #Try again with the first alignment of source
                        catfile=open(tmpsa,'w')
                        for f in initstran:
                            catfile.write(f)
                        catfile.close()
                        if(verbose==1):
                            print("next round will start again with the initial alignment score= %s" % (gfithfit[best,0]))
                            mfile.write("next round will start again with the initial alignment score= %s \n" % (gfithfit[best,0]))

            nbri=nbri+1

#print results of best
#mfile=open(output, 'a')
mfile.write("Best solution (Source_tran.sdf):  \n")
mfile.write("gfit= %s cfit= %s hfit= %s gfit+hfit= %s \n" % (gfithfit[best,1],gfithfit[best,3],gfithfit[best,2],gfithfit[best,0]))
mfile.close()
if(verbose==1):
    print("")
    print("Best solution (Source_tran.sdf):")
    print("gfit= %s cfit= %s hfit= %s gfit+hfit= %s" % (gfithfit[best,1],gfithfit[best,3],gfithfit[best,2],gfithfit[best,0]))
if(isexist("Source_tran.sdf")):
        os.remove("Source_tran.sdf")
catfile=open("Source_tran.sdf", 'a')
for f in stran:
    catfile.write(f)
catfile.close()

if(flagcleanfiles==1):
    if(isexist("tmpt-dots.pdb")):
        os.remove("tmpt-dots.pdb")
    if(isexist("tmpsf.sdf")):
        os.remove("tmpsf.sdf")
    if(isexist("tmptf.sdf")):
        os.remove("tmptf.sdf")
    if(isexist("tmpsa.sdf")):
        os.remove("tmpsa.sdf")
    if(isexist("Source_tran-conf.sdf")):
        os.remove("Source_tran-conf.sdf")
    if(isexist("make_fgts.sdf")):
        os.remove("make_fgts.sdf")
    ri=1
    while (ri <= nbfgt):
        name="make_fgts_" + str(ri) + ".sdf"
        named="Source_" + str(ri) + "-dots.pdb"
        namedc="Sourceconf_" + str(ri) + "-dots.pdb"
        if(isexist(name)):
            os.remove(name)
        if(isexist(named)):
            os.remove(named)
        if(isexist(namedc)):
           os.remove(namedc)
        ri=ri+1

    if(isexist("combination.sdf")):
        os.remove("combination.sdf")
    if(isexist("combination-min.sdf")):
        os.remove("combination-min.sdf")
    if(isexist("slog")):
        os.remove("slog")
    if(isexist("tran.txt")):
        os.remove("tran.txt")

if(verbose==1):
    print("Done (see slogflex)")

##################################
