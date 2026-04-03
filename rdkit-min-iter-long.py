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

#m = Chem.SDMolSupplier(molecule,removeHs=True)
m = Chem.SDMolSupplier(molecule,removeHs=False)
w = Chem.SDWriter(output)

for mol in m:
    Chem.AssignAtomChiralTagsFromStructure(mol)
    
    #with UFF removeHs=True is ok
    #ff = AllChem.UFFGetMoleculeForceField(mol, confId=-1)
   
    #with MMFF94 
    mmff_props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant='mmff')
    ff = AllChem.MMFFGetMoleculeForceField(mol, mmff_props, confId=-1, ignoreInterfragInteractions=True)

    einit = ff.CalcEnergy();
    #print('Einit %s' % einit)
    count=1
    diff=100
    eold = einit
    #while (diff > 0.001) :
    #while count <= 500 and count2 < 10 :
    #while count <= 500 :
    while count <= 2000 :
        ff.Minimize(1)
        e = ff.CalcEnergy();
        diff = eold - e
        eold = e
        #print('iter %s = %s (diff= %s)' % (count, e, diff))
        count +=1
    #e = ff.CalcEnergy();
    print('Einit %s to Eend %s' % (einit, eold))
    
    #mol=Chem.AddHs(mol,addCoords=True)
    mol=Chem.AddHs(mol,addCoords=False)

    w.write(mol,confId=-1)

w.flush()
w.close()

#######################################
