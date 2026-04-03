#!/usr/bin/python3.7

#execute ./<>.py sdf_file sdf_outputname  

import os, sys
from rdkit import Chem
from rdkit.Chem import AllChem

# sys.argv[0] is the name of the program itself
molecule=sys.argv[1]
output=sys.argv[2]

#######################################
# MAIN program

#if set True: problems occur with kekulization:
m = Chem.SDMolSupplier(molecule,removeHs=False)
w = Chem.SDWriter(output)

for mol in m:
    #if de-commented: problem occur when combination.sdf is weird
    Chem.AssignAtomChiralTagsFromStructure(mol)

    #with UFF removeHs=True is ok
    #ff = AllChem.UFFGetMoleculeForceField(mol, confId=-1)
   
    #with MMFF94 
    mmff_props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant='mmff')
    ff = AllChem.MMFFGetMoleculeForceField(mol, mmff_props, confId=-1, ignoreInterfragInteractions=True)

    einit = ff.CalcEnergy();
    #print('Einit %s' % einit)
    count=1
    count2=0
    diff=100
    eold = einit
    #while count <= 1000 :
    while count <= 500 :
    #while count <= 800 :
        ff.Minimize(1)
        e = ff.CalcEnergy();
        diff = eold - e
        eold = e
        #print('iter %s = %s (diff= %s)' % (count, e, diff))
        count +=1
        if diff < 1 :
            count2 +=1
        else :
            count2 = 0

    #e = ff.CalcEnergy();
    #print('Einit %s to Eend %s' % (einit, eold))
    
    mol=Chem.AddHs(mol,addCoords=False)
    w.write(mol,confId=-1)

w.flush()
w.close()

#######################################
