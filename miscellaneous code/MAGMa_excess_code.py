#!/usr/bin/env python3
"""
Author: Anna-Marie Seelen
Studentnumber:1008970
Description: alignment of protein sequences to a reference protein sequence
and returning stats
Usage: python3 *name_of_script* *name_input_file* [*output_file*] [gap_penalty]
    name_of_script:
    name_input_file:
    output_file: (optional)
    gap_penalty: default is 8 (optional)
"""

# import statements
import json
import sqlite3
from sys import argv
import re
from pathlib import Path
import os
import subprocess
import sys
from decimal import *
import ast
from matchms import Spectrum
from matchms.filtering import add_losses
from rdkit import Chem
from rdkit.Chem import rdFMCS

def convert_json_to_sqlite(path_to_json_file,dir_database):
    connection = sqlite3.connect('{0}'.format(dir_database))
    cursor = connection.cursor()
    try:
        cursor.execute("""CREATE TABLE molecules (id TEXT PRIMARY KEY,
                                             mim INTEGER NOT NULL,
                                             charge INTEGER NOT NULL,
                                             natoms INTEGER NOT NULL,
                                             molblock TEXT,
                                             inchikey TEXT,
                                             smiles TEXT,
                                             molform TEXT,
                                             name TEXT,
                                             reference TEXT,
                                             logp INT)""")
        connection.commit()
        print("HMDB_MAGMa.db created")
    except:
        print ("HMDB_MAGMa.db already exists (or error creating it)")
        exit()

    HMDB=json.load(open(path_to_json_file))
    columns = ['spectrum_id', 'peaks_json', 'Ion_Mode', 'Compound_Name', 'Precursor_MZ', 'ExactMass', 'Charge', 'Smiles']
    #problem: the peak formatting is not in mgf style...
    for row in HMDB:
        #keys = tuple(row[c] for c in columns)
        #cursor.execute('insert into Spectrum values(?,?,?,?,?,?,?,?)', keys)
        cursor.execute('''INSERT INTO molecules (id, mim, charge, natoms, molblock, inchikey,
                                     smiles,molform,name,reference,logp) VALUES (?,?,?,?,?,?,?,?,?,?,?)''', (
            hmdb_id,
            int(mim * 1e6),
            charge,
            int(natoms),
            unicode(molblock),
            unicode(inchikey),
            unicode(smiles),
            unicode(molform),
            unicode(molname, 'utf-8', 'xmlcharrefreplace'),
            unicode(hmdb_id),
            int(logp * 10)))
        memstore[inchikey] = (hmdb_id, hmdb_id, ionized)
    conn.commit()
        #print(f'{row["spectrum_id"]} data inserted Succesfully')

    cursor.execute(
        'CREATE INDEX idx_cover ON molecules (charge,mim,natoms,reference,molform,inchikey,smiles,name,molblock,logp)')
    connection.commit()
    return None

def select_features(list_of_features,list_with_fragments_and_smiles):
    for feature in list_of_features:
        if re.search(r'loss', feature[0]) != None:
            # first calculate all the losses in de spectrum based on the parent mass
            spectrum = add_losses(spectrum)
            for loss in range(len(spectrum.losses.mz)):
                # get the losses with 2 decimal points and round down the number
                rounded_loss = Decimal(spectrum.losses.mz[loss]).quantize(Decimal('.01'),
                                                                          rounding=ROUND_DOWN)
                # if a loss in the spectrum contains the loss in the feature then
                # the document contains the loss!
                if float(rounded_loss) == float(re.search(r'\_(.*\..{2}).*', feature[0]).group(1)):
                    num_of_ass_doc_with_feature += 1
                    documents_that_contain_feature.append(int(document[0]))
    feature_number = re.search(r'\_(.*)', object[0]).group(1)
    rounded_feature = Decimal(feature_number).quantize(Decimal('.01'), rounding=ROUND_DOWN)
    list_of_features.append('%.2f' % float(rounded_feature))
    #if molid is not None, than you can look for mol id in the fragments table
    return None
#you need the spectrum file with the motif and the identifier
#you need the file from pdf script where you see the fragments for each motif

def search_for_smiles(list_of_features,list_with_fragments_and_smiles):
    # for every feature in the motif
    for feature in list_of_features:
        # if the feature in the motif is a loss
        if re.search(r'loss', feature) != None:
            # make a list of losses from the list_with_fragments_and_smiles from MAGMa
            list_of_losses=make_list_of_losses(list_with_fragments_and_smiles)
            # and for every loss in this list of losses, try to find a match to the loss of the feature
            for index,rounded_loss in enumerate(list_of_losses):
                # if you did find a match up to 2 decimals
                if float(rounded_loss) == float(re.search(r'\_(.*\..{2}).*', feature).group(1)):
                    # then retrieve the corresponding fragment string and the parent string beloning to the loss
                    precusor_mz, precusor_smiles = list_with_fragments_and_smiles[0]
                    fragment_mz, fragment_smiles= list_with_fragments_and_smiles[index]
                    smiles_neutral_loss=make_regex_for_loss(precusor_smiles, fragment_smiles)
                    # check if the smiles has the same molecular mass as the loss reported of the feature
                    mol_weight_from_smiles=MolWt(Chem.MolFromSmiles(f'{smiles_neutral_loss}'))
                    rounded_mol_weight_from_smiles=Decimal(mol_weight_from_smiles).quantize(Decimal('.1'),
                                                                          rounding=ROUND_DOWN)
                    if rounded_mol_weight_from_smiles==float(re.search(r'\_(.*\..{1}).*', feature).group(1)):
                        print(smiles_neutral_loss)
        # if feature is not a loss
        else:
            for fragment in enumerate(list_with_fragments_and_smiles):
                fragment_mz, fragment_smiles = fragment
                rounded_fragment = Decimal(fragment_mz).quantize(Decimal('.01'),
                                                                          rounding=ROUND_DOWN)
                if float(rounded_fragment)==float(re.search(r'\_(.*\..{2}).*', feature).group(1)):
                    print(fragment_smiles)
    return None

                    # then we can look for the smiles! with the parent

    #             feature_loss = re.search(r'\_(.*)', feature).group(1)
    #             loss_in_feature_list= Decimal(feature_loss).quantize(Decimal('.01'),
    #                                                                   rounding=ROUND_DOWN)
    #             if calculated_loss==loss_in_feature_list:
    #
    #         spectrum = add_losses(spectrum)
    #     for loss in range(len(spectrum.losses.mz)):
    #         # get the losses with 2 decimal points and round down the number
    #         rounded_loss = Decimal(spectrum.losses.mz[loss]).quantize(Decimal('.01'),
    #                                                                   rounding=ROUND_DOWN)
    #         # if a loss in the spectrum contains the loss in the feature then
    #         # the document contains the loss!
    #         feature_number = re.search(r'\_(.*)', feature).group(1)
    #         rounded_feature = Decimal(feature_number).quantize(Decimal('.01'), rounding=ROUND_DOWN)
    #         #rounded_feature='%.2f' % float(rounded_feature)
    #         if float(rounded_loss) == float(re.search(r'\_(.*\..{2}).*', feature[0]).group(1)):
    #             num_of_ass_doc_with_feature += 1
    #             documents_that_contain_feature.append(int(document[0]))
    # # if the feature we are looking at is a fragment
    # else:
    #     for fragment in range(len(spectrum.peaks.mz)):
    #         # get the fragment of the spectrum with 2 decimal points and round down the number
    #         rounded_fragment = Decimal(spectrum.peaks.mz[fragment]).quantize(Decimal('.01'),
    #                                                                          rounding=ROUND_DOWN)
    #         # if a fragment in the spectrum is the same as a fragment in the feature
    #         # then the document contains the fragment!
    #         if float(rounded_fragment) == float(re.search(r'\_(.*\..{2}).*', feature[0]).group(1)):
    #             print("hi")

    #alle haakjes weg
    #only the letters should match, if there is a ( I don't care) regex expression
    #string_1='NC(CC(=O)O)C(=O)O'
    #string_2='N[\(|\)]*C[\(|\)]*C[\(|\)]*C[\(|\)]*\=[\(|\)]*O[\(|\)]*\)[\(|\)]*O[\(|\)]*'
    #mols = [Chem.MolFromSmiles(string_2), Chem.MolFromSmiles(string_1)]
    #MCS_in_smart_string=rdFMCS.FindMCS(mols, ringMatchesRingOnly=True).smartsString
    #print(Chem.MolToSmiles(Chem.MolFromSmarts(MCS_in_smart_string)))
    #print(re.search(fr'(.*)({0})(.*)'.format(string_2), string_1).group(1))
    #string_2="NCCC(=O)O"
    #print(re.search(r'(.*)(N[\(|\)]*C[\(|\)]*C[\(|\)]*C[\(|\)]*\=[\(|\)]*O[\(|\)]*\)[\(|\)]*O[\(|\)]*)(.*)', string_1).group(3))