# Script to convert chemicals provided in the molecular representation wihin a SDF file to InChIKeys
# Workflow is based on the ChEMBL structure curation pipeline using RDKit (see Readme for details)


import rdkit.Chem as Chem
import chembl_structure_pipeline
import csv
import os
import string
import random
import itertools
import pandas as pd
import numpy as np

# This exmaple script is written to deal with SDF files as provided by the commercial vendor Selleck chemicals.
#        Adaptions will be neccessary if a SDF from a different source is used. 

# Step 1: 
# read the molecules from the SDF file using a context manager 

# additional notes: 
# for each element of the SDF file check whether a molecule is present and whether is has at least one atom
# - for each molecule we read the annotation file which is saved in a dictionary (named list)
# - we save the molecules in a dictionary with the vnedor's catalog number as keys 
# - for some molecules the entries cannot ba read because the IDs contain non Unicode characters and we assign dummy names to them

my_molecule_dict = {}
count=0

with Chem.SDMolSupplier('./my_sdf.SDF') as sdf:
        for mol in sdf:
            if mol is None: continue
            if mol.GetNumAtoms() > 0: 
                count +=1
                try:
                    molecule_dict = mol.GetPropsAsDict()
                    mol_keys = list(molecule_dict.keys())
                    mol_keys_cn = ["Catalog Number" in i for i in mol_keys]
                    mol_key = list(itertools.compress(mol_keys, mol_keys_cn))[0]
                    mol_cn = molecule_dict[mol_key]
                    my_molecule_dict[mol_cn] = mol
                except (UnicodeDecodeError,IndexError):
                    mol_cn = "cn_not_read" + str(count)
                    my_molecule_dict[mol_cn] = mol
                    
# extract the molecules from the dictionary as a list
molecule_list = list(my_molecule_dict.values())


# check if and how many molecules were not read in 
SDF_full = Chem.SDMolSupplier('./my_sdf.SDF')
print(str(len(SDF_full)-len(molecule_list))+" molecules not read")
                           


# Step 2: 
# standardize the molecules (i.e. strip salts, remove charges)

from chembl_structure_pipeline import standardizer
molecule_list = list(map(standardizer.standardize_mol, molecule_list))


# Step 3: 
# get parent molecule 
# note: this resturns a tuple for which the first entry contains the molecules (standardizer is imported from chembl_structure_pipeline)

molecule_list = list(map(standardizer.get_parent_mol, molecule_list))
molecule_list = [x[0] for x in molecule_list]


# Step 4:
# convert molecules to InChI and InChIKey

inchi_from_mol = list(map(Chem.rdinchi.MolToInchi, molecule_list))
inchi_from_mol = [x[0] for x in inchi_from_mol]
inchikey_from_inchi = list(map(Chem.rdinchi.InchiToInchiKey, inchi_from_mol))

# Step 4:
# save the InChIKeys with the catalog numbers as csv 
df_out = pd.DataFrame({'InChIKey':inchikey_from_inchi})
df_out['catalog_number'] = list(my_molecule_dict.keys())


df_out.to_csv("./InChIKey_standardized_fromSDF.csv")
