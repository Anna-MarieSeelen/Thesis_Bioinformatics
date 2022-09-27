#!/usr/bin/env python3
"""
Author: Anna-Marie Seelen
Studentnumber:1008970
Description: This script runs MAGMa on one mgf formatted spectrum file and looks for smiles annotations for the features
of the corresponding mass2motif for which the spectrum file was made using MassQL
Usage: python3 *name_of_script* *path_to_structures_database* *path_to_spectrum_file* *path_to_store_results_db*
*path_to_txt_file_with_motif_and_frag*

    path_to_structures_database: path where the structure database is which was downloaded from HMDB and formatted with
    process_hmdb.py from MAGMa
    path_to_spectrum_file: path where the spectrum file is stored with a spectrum from the HMDB that contains a motif
    which is determined with MassQL.py script.
    path_to_store_results_db: path where the output database from MAGMa with the annotations for the spectrum should be
    stored
    path_to_txt_file_with_motif_and_frag: path to the output file from make_pdf_with_smiles.py in which each feature of
    each mass2motif is indicated.
"""

# import statements
import sqlite3
from sys import argv
import re
from pathlib import Path
import os
import subprocess
from decimal import *
import ast
from rdkit import Chem
from rdkit.Chem.Descriptors import MolWt

#functions
def initialize_db_to_save_results(path_to_store_results_db,path_to_spectrum_file):
    """
    Initializes a sqlite database using the MAGMa CL init_db function to store annotation results from MAGMa for spectrum

    :param path_to_store_results_db: str, path where the output database from MAGMa with the annotations for the
    spectrum should be stored.
    :param path_to_spectrum_file: str, path where the spectrum file is stored with a spectrum from the HMDB that
    contains a motif which was determined with MassQL.py script.
    :return: str containing the realpath to the results database from MAGMa
    """
    #generate a good database name with the spectrum id and motif in it.
    spectrum_name=re.search(r'(spectrum.*)(HMDB.*)(.txt)', path_to_spectrum_file).group(2)
    #realpath of the database
    file_path_out = Path(r"{0}/MAGMa_db_{1}.sqlite".format(path_to_store_results_db, spectrum_name))
    #if you try to annotate something in an existing database it will go wrong, so if the database exists remove it.
    if os.path.exists(file_path_out):
        os.remove(file_path_out)
    cmd = 'magma init_db {0}'\
            .format(file_path_out)
    try:
        e = subprocess.check_call(cmd, shell=True, stdout = subprocess.DEVNULL, stderr = subprocess.STDOUT)
    # For some reason the text that MAGMa returns on the command line is seen as an error so include this line to keep
    # script from stopping.
    except subprocess.CalledProcessError:
        pass
    return file_path_out

def add_spectrum_to_be_annotated_into_db(path_to_results_db_file, path_to_spectrum_file, abs_intensity_thres=1000, mz_precision_ppm=80, mz_precision_abs=0.01, spectrum_file_type='mgf', ionisation=1):
    """
    Adds the spectrum to be annotated in the sqlite results database using the MAGMa CL read_ms_data function.

    :param path_to_results_db_file: str, path to the results database in which the spectrum and annotations will be
    stored
    :param path_to_spectrum_file: str, path where the spectrum file is stored with a spectrum from the HMDB that
    contains a motif which is determined with MassQL.py script.
    :param abs_intensity_thres: int, absolute intensity threshold for storing peaks in database (default: 1000)
    :param mz_precision_ppm: int, maximum relative m/z error (ppm) used to determine if structure in structure database
    is good fit (default: 80)
    :param mz_precision_abs: float, maximum absolute m/z error (Da) used to determine if structure in structure database
    is good fit (default: 0.01)
    :param spectrum_file_type: str, indicates how the spectrum_file is formatted (default: mgf)
    :param ionisation: {1, -1}, indicates the ionisation mode, -1 is negative, 1 is positive (default: 1)
    :return: None
    """
    cmd = 'magma read_ms_data -f {0} -i {1} -a {2} -p {3} -q {4} {5} {6}'.format(spectrum_file_type,ionisation,abs_intensity_thres, mz_precision_ppm, mz_precision_abs, path_to_spectrum_file, path_to_results_db_file)
    e = subprocess.check_call(cmd, shell=True,stdout = subprocess.DEVNULL, stderr = subprocess.STDOUT)
    return None

def annotate_spectrum_with_MAGMa(path_to_structures_database,path_to_results_db_file, max_num_break_bonds=10, structure_db="hmdb", ncpus=1):
    """
    Finds matches for spectrum in structure database and annotates generated fragments using MAGMa CL annotate function.

    :param path_to_structures_database: str, path where the structure database is which was downloaded from HMDB and
    formatted with process_hmdb.py from MAGMa
    :param path_to_results_db_file: str, path to the results database in which the spectrum and annotations will be
    stored
    :param max_num_break_bonds: int, maximum number of bond breaks to generate substructures (default: 10)
    :param structure_db: {pubchem,kegg,hmdb}, structure database from which the matches should be retrieved
    (default: hmdb)
    :param ncpus: int, number of parallel cpus to use for annotation (default: 1)
    :return: None
    """
    cmd = "magma annotate -b {0} -s {1} -o {2} -n {3} {4}".format(max_num_break_bonds,structure_db,path_to_structures_database,ncpus, path_to_results_db_file)
    e = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
    return None

def get_molid_of_matched_compound(path_to_results_db_file, identifier_spectrum):
    """
    Retrieves the molid for the compound with the same HMDB identifier as the HMDB identifier of the spectrum file.

    :param path_to_results_db_file: str, path to the results database in which the spectrum and annotations are stored
    :param identifier_spectrum: str, the HMDB identifier of the spectrum
    :return: None or int, the molid for the compound with the same HMDB identifier as the HMDB identifier of the
    spectrum file.

    Basically you only want the MAGMa annotation for the compound which is the same compound as the spectrum, and since
    we are using the same database for the structures and for MassQL we can look up the right compound with the HMDB
    identifier.
    """
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
        #if the identifier of the spectrum matches an identifier in the results database, fetch that molid
        if identifier_spectrum==str(re.search(r'.*\((HMDB.*)\).*', identifier).group(1)):
            conn = sqlite3.connect(path_to_results_db_file)
            sqlite_command = \
                f"""SELECT molid FROM molecules WHERE name = '{identifier}'"""
            cur = conn.cursor()
            cur.execute(sqlite_command)
            molid=cur.fetchall()
            molid=[i[0] for i in molid][0]
            return molid
        #if the identifier of the spectrum does not match an identifier in the results database, MAGMa did not find the
        # compound the spectrum actually belongs to, so return None, we cannot use the fragment annotations.
        else:
            return None

def fetch_fragments_and_annotations_for_molid(molid,path_to_results_db_file):
    """
    Makes list of tuples containing the fragment m/z and smiles for each fragment of indicated molid

    :param molid: int, the molid for the compound with the same HMDB identifier as the HMDB identifier of the spectrum
    file.
    :param path_to_results_db_file: str, path to the results database in which the spectrum and annotations are stored.
    :return: list of tuples containing the fragment m/z and smiles for each fragment of indicated molid
    """
    conn = sqlite3.connect(path_to_results_db_file)
    sqlite_command = \
        f"""SELECT mz,smiles FROM fragments WHERE molid = '{molid}'"""
    cur = conn.cursor()
    cur.execute(sqlite_command)
    list_with_fragments_and_smiles=cur.fetchall()
    return list_with_fragments_and_smiles

def parse_line_with_motifs_and_querries(line):
    """
    Takes a line in a tab separated file and separates it into motif, feature list and massql query
    :param line: str, line in a tab separated file containing the motif, list of features and the massql query
    :return: returns the motif, feature_list and query in a tuple
    """
    motif=re.search(r'(.*)    (.*)    (.*)', line).group(1)
    feature_list=re.search(r'(.*)    (.*)    (.*)', line).group(2)
    query=re.search(r'(.*)    (.*)    (.*)', line).group(3)
    return motif,feature_list,query

def look_for_features(txt_file_with_motif_and_frag, current_motif):
    """
    Returns the list of features that make up the motif that are according to MassQL present in the spectrum file.

    :param txt_file_with_motif_and_frag: str, path to the output file from make_pdf_with_smiles.py in which each feature
    of each mass2motif is indicated.
    :param current_motif: str, motif in the name of the spectrum file, so the motif, which is according to MassQL,
    present in the spectrum.
    :return: list of feature beloning to the mass2motif which is according to MassQL present in the spectrum.
    """
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
    """

    :param list_with_fragments_and_smiles:
    :return:
    """
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

def search_for_smiles(list_of_features,list_with_fragments_and_smiles):
    """

    :param list_of_features:
    :param list_with_fragments_and_smiles:
    :return:
    """
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
    #TODO: fix output here
    return None

def write_output_to_file():
    # each time we will write a annotation for a similar fragment in the output file
    # tab separated file with motif feature annotation_1 annotation_2 annotation_3
    # TODO: for one feature you might find multiple annotations
    # TODO: kijken naar welke dingen je wil returnen en hoe de output van dit hele script eruit moet zien
    return None

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
    #step 1: initialize database to save results from MAGMa
    path_to_results_db_file=initialize_db_to_save_results(path_to_store_results_db, path_to_spectrum_file)
    # step 2: add spectrum to be annotated into the results database
    add_spectrum_to_be_annotated_into_db(path_to_results_db_file, path_to_spectrum_file, abs_intensity_thres=1000,
                                         mz_precision_ppm=80, mz_precision_abs=0.01, spectrum_file_type='mgf',
                                         ionisation=1)
    #step 3: annotate spectrum with MAGMa and store output in results database
    annotate_spectrum_with_MAGMa(path_to_structures_database, path_to_results_db_file, max_num_break_bonds=10, structure_db="hmdb",
                      ncpus=1)
    # step 4: get identifiers of matches
    molid=get_molid_of_matched_compound(path_to_results_db_file, identifier)
    #step 5: search for the features beloning to motif
    list_of_features=look_for_features(path_to_txt_file_with_motif_and_frag, current_motif)
    print(list_of_features)
    # step 6: get the fragments and smiles of the match
    if molid is not None:
        list_with_fragments_and_smiles=fetch_fragments_and_annotations_for_molid(molid, path_to_results_db_file)
        search_for_smiles(list_of_features,list_with_fragments_and_smiles)
        #search_for_smiles(list_of_features, list_with_fragments_and_smiles)
    # step 7:

if __name__ == "__main__":
    main()