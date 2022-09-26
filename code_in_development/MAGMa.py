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
from decimal import *
import ast
from matchms.filtering import add_losses
from rdkit import Chem
from rdkit.Chem import rdFMCS
from rdkit.Chem.Descriptors import MolWt

#functions
def initialize_db_to_save_results(path_to_store_results_db,path_to_spectrum_file):
    #generate a good database name with the spectrum id and motif in it.
    spectrum_name=re.search(r'(.*)(spectrum.*)(.txt)', path_to_spectrum_file).group(2)
    file_path_out = Path(r"{0}/MAGMa_db_{1}.sqlite".format(path_to_store_results_db, spectrum_name))
    #if you try to annotate something in an existing database it will go wrong, so if the database exists remove it.
    if os.path.exists(file_path_out):
        os.remove(file_path_out)
    cmd = 'magma init_db {0}'\
            .format(file_path_out)
    try:
        e = subprocess.check_call(cmd, shell=True, stdout = subprocess.DEVNULL, stderr = subprocess.STDOUT)
    except subprocess.CalledProcessError:
        pass
    return file_path_out

def add_spectrum_to_be_annotated_into_db(path_to_results_db_file, path_to_spectrum_file, abs_intensity_thres=1000, mz_precision_ppm=80, mz_precision_abs=0.01, spectrum_file_type='mgf', ionisation=1):
    # a = Absolute intensity threshold for storing peaks in database (default: 1000)
    # i = Ionisation mode (default: 1)
    # p= Maximum relative m/z error (ppm) (default: 5)
    # q= Maximum absolute m/z error (Da) (default: 0.001)
    cmd = 'magma read_ms_data -f {0} -i {1} -a {2} -p {3} -q {4} {5} {6}'.format(spectrum_file_type,ionisation,abs_intensity_thres, mz_precision_ppm, mz_precision_abs, path_to_spectrum_file, path_to_results_db_file)
    e = subprocess.check_call(cmd, shell=True,stdout = subprocess.DEVNULL, stderr = subprocess.STDOUT)
    #print("EXIT STATUS AND TYPE", e, type(e))
    return None

def annotate_spectrum(path_to_structures_database,path_to_results_db_file, max_num_break_bonds=10, structure_db="hmdb", ncpus=1):
    # b = Maximum number of bond breaks to generate substructures (default: 3)
    # s= --structure_database {pubchem,kegg,hmdb} Retrieve molecules from structure database (default: )
    # n = Number of parallel cpus to use for annotation (default: 1)
    cmd = "magma annotate -b {0} -s {1} -o {2} -n {3} {4}".format(max_num_break_bonds,structure_db,path_to_structures_database,ncpus, path_to_results_db_file)
    e = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
    # e=e.decode('utf-8')
    # for i in e.split('\n'):
    #     print(i)
    # if you want to do something with the output
    return None

def get_molid_of_matched_compound(path_to_results_db_file, identifier_spectrum):
    conn = sqlite3.connect(path_to_results_db_file)
    sqlite_command = \
        f"""SELECT name FROM molecules"""
    cur = conn.cursor()
    cur.execute(sqlite_command)
    list_of_identifiers_of_matches = cur.fetchall()
    list_of_identifiers_of_matches = list(map(''.join,list_of_identifiers_of_matches))
    print(list_of_identifiers_of_matches)
    for identifier in list_of_identifiers_of_matches:
        print(str(re.search(r'.*\((HMDB.*)\).*', identifier).group(1)))
        if identifier_spectrum==str(re.search(r'.*\((HMDB.*)\).*', identifier).group(1)):
            conn = sqlite3.connect(path_to_results_db_file)
            sqlite_command = \
                f"""SELECT molid FROM molecules WHERE name = '{identifier}'"""
            cur = conn.cursor()
            cur.execute(sqlite_command)
            molid=cur.fetchall()
            molid=[i[0] for i in molid][0]
            return molid
        else:
            pass
            return None

def make_spectrum_object(molid,path_to_results_db_file):
    conn = sqlite3.connect(path_to_results_db_file)
    sqlite_command = \
        f"""SELECT mz,smiles FROM fragments WHERE molid = '{molid}'"""
    cur = conn.cursor()
    cur.execute(sqlite_command)
    list_with_fragments_and_smiles=cur.fetchall()
    return list_with_fragments_and_smiles

def parse_line_with_motifs_and_querries(line):
    motif=re.search(r'(.*)    (.*)    (.*)', line).group(1)
    feature_list=re.search(r'(.*)    (.*)    (.*)', line).group(2)
    query=re.search(r'(.*)    (.*)    (.*)', line).group(3)
    return motif,feature_list,query

def look_for_features(txt_file_with_motif_and_frag, current_motif):
    lines = (open(txt_file_with_motif_and_frag))
    print("opened file")
    for line in lines:
        line = line.strip()
        line = line.replace('\n', '')
        motif, features_list, query = parse_line_with_motifs_and_querries(line)
        features_list_1 = ast.literal_eval(features_list)
        if motif == current_motif:
            list_of_features=[]
            for object in features_list_1:
                list_of_features.append(object[0])
            return(list_of_features)

def make_list_of_losses(list_with_fragments_and_smiles):
    # make a list of losses if re.search is a loss
    list_with_losses=[]
    for index,fragment in enumerate(list_with_fragments_and_smiles):
        if index !=0:
            precusor_mz, precusor_smiles=list_with_fragments_and_smiles[0]
            fragment_mz, fragment_smiles=fragment
            calculated_loss=Decimal(precusor_mz-fragment_mz).quantize(Decimal('.01'),
                                                                      rounding=ROUND_DOWN)
            list_with_losses.append(calculated_loss)
        else:
            precusor_mz, precusor_smiles = fragment
            list_with_losses.append(precusor_mz)
    return(list_with_losses)

def search_for_smiles(list_of_features,list_with_fragments_and_smiles):
    #TODO: kijken naar welke dingen je wil returnen en hoe de output van dit hele script eruit moet zien
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

def make_regex_for_loss(parent_string,fragment_string):
    #TODO: werkt dit voor alle moleculen? Aan justin vragen...
    new_fragment_string=re.escape(fragment_string)
    list_fragment_string = list(new_fragment_string)
    for index,letter in enumerate(list_fragment_string):
        if index==0:
            if letter=="\\":
                letter_with_slash="[\(|\)]*"+letter
                list_fragment_string[index]=letter_with_slash
            else:
                letter_with_slash="[\(|\)]*"+letter+"[\(|\)]*"
                list_fragment_string[index]=letter_with_slash
        elif letter!="\\":
            if index!=0:
                letter_with_slash=letter+"[\(|\)]*"
                list_fragment_string[index]=letter_with_slash
    re_search_string="".join(list_fragment_string)
    if re.search(rf'(.+)({re_search_string})',
                    parent_string) != None:
        smiles_neutral_loss = re.search(rf'(.+)({re_search_string})',
                    parent_string).group(1)
        return smiles_neutral_loss
    elif re.search(rf'({re_search_string})(.+)',
                    parent_string) != None:
        smiles_neutral_loss = re.search(rf'({re_search_string})(.+)',
                                        parent_string).group(2)
        return smiles_neutral_loss
    else:
        return "smiles not found"

def main():
    #main function of the script
    #step 0: parse input
    path_to_structures_database=argv[1]
    path_to_spectrum_file=argv[2]
    path_to_store_results_db=argv[3]
    path_to_txt_file_with_motif_and_frag=argv[4]
    # this script is made to return the smiles corresponding to the features of a motif of one spectrum_file_from_massql
    # after parsing the name of the spectrum file you will get this:
    identifier="HMDB0000191"
    current_motif="gnps_motif_38.m2m"
    #step 1: initialize database to save results
    path_to_results_db_file=initialize_db_to_save_results(path_to_store_results_db, path_to_spectrum_file)
    # step 2: add spectrum to be annotated into the results database
    add_spectrum_to_be_annotated_into_db(path_to_results_db_file, path_to_spectrum_file, abs_intensity_thres=1000,
                                         mz_precision_ppm=80, mz_precision_abs=0.01, spectrum_file_type='mgf',
                                         ionisation=1)
    #step 3: annotate spectrum and store output in results database
    annotate_spectrum(path_to_structures_database, path_to_results_db_file, max_num_break_bonds=10, structure_db="hmdb",
                      ncpus=1)
    # step 4: get identifiers of matches
    molid=get_molid_of_matched_compound(path_to_results_db_file, identifier)
    #step 5: search for the features beloning to motif
    list_of_features=look_for_features(path_to_txt_file_with_motif_and_frag, current_motif)
    print(list_of_features)
    # step 6: get the fragments and smiles of the match
    if molid is not None:
        list_with_fragments_and_smiles=make_spectrum_object(molid, path_to_results_db_file)
        search_for_smiles(list_of_features,list_with_fragments_and_smiles)
        #search_for_smiles(list_of_features, list_with_fragments_and_smiles)
    # step 7:

if __name__ == "__main__":
    main()