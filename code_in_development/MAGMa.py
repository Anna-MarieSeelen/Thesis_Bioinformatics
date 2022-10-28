#!/usr/bin/env MAGMa
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
import time
import shutil
import pandas as pd
import numpy as np

#functions

def construct_path_to_db(path_to_spectrum_file: str, path_to_store_results_db: str) -> tuple:
    # generate a good database name with the spectrum id and motif in it.
    print(path_to_spectrum_file)
    path, filename = os.path.split(path_to_spectrum_file)
    print(filename)
    identifier = re.search(r'(spectrum_)(.*motif.*)_(HMDB.*)(.txt)', filename).group(3)
    motif = re.search(r'(spectrum_)(.*motif.*)_(HMDB.*)(.txt)', filename).group(2)
    # realpath of the database
    file_path_out = Path(fr"{path_to_store_results_db}/MAGMa_db_{identifier}.sqlite")
    # if you try to annotate something in an existing database it will go wrong, so if the database exists do not
    # annotate it again.
    if os.path.exists(file_path_out):
        #assert False, "The results db for this spectrum already exists, remove it"
        return "exists", file_path_out, identifier, motif
    else:
        return "new", file_path_out, identifier, motif


def initialize_db_to_save_results(path_to_results_db) -> None:
    """
    Initializes a sqlite database using the MAGMa CL init_db function to store annotation results from MAGMa for spectrum

    :param path_to_store_results_db: str, path where the output database from MAGMa with the annotations for the
    spectrum should be stored.
    :param path_to_spectrum_file: str, path where the spectrum file is stored with a spectrum from the HMDB that
    contains a motif which was determined with MassQL.py script.
    :return: str containing the realpath to the results database from MAGMa
    """

    cmd = f'magma init_db {path_to_results_db}'
    try:
        e = subprocess.check_call(cmd, shell=True, stdout = subprocess.DEVNULL, stderr = subprocess.STDOUT)
    # For some reason the text that MAGMa returns on the command line is seen as an error so include this line to keep
    # script from stopping.
    except subprocess.CalledProcessError:
        pass
    return None

def add_spectrum_into_db(path_to_results_db_file: str, path_to_spectrum_file: str,
                                             abs_intensity_thres=1000, mz_precision_ppm=80, mz_precision_abs=0.01,
                                             spectrum_file_type='mgf', ionisation=1) -> None:
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
    cmd = f'magma read_ms_data -f {spectrum_file_type} -i {ionisation} -a {abs_intensity_thres} -p {mz_precision_ppm} -q {mz_precision_abs} {path_to_spectrum_file} {path_to_results_db_file}'
    e = subprocess.check_call(cmd, shell=True,stdout = subprocess.DEVNULL, stderr = subprocess.STDOUT)
    return None

#TODO: add_structures and remove structure database etc.

magma add_structures -t smiles 'C(C(C(=O)O)N)C(O)=O' results.sqlite

def annotate_spectrum_with_MAGMa(path_to_structures_database: str, path_to_results_db_file: str,
                                     max_num_break_bonds=10, structure_db="hmdb", ncpus=1) -> None:
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
    cmd = f"magma annotate -b {max_num_break_bonds} -s {structure_db} -o {path_to_structures_database} -n {ncpus} {path_to_results_db_file}"
    e = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
    return None

#TODO: deze functie kan verwijderd worden, want je hebt nu maar 1 molid.

def get_molid_of_matched_compound(path_to_results_db_file: str, identifier_spectrum: str) -> int:
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
    print(path_to_results_db_file)
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
            print(molid)
            return molid
        #if the identifier of the spectrum does not match an identifier in the results database, MAGMa did not find the
        # compound the spectrum actually belongs to, so return None, we cannot use the fragment annotations.
    print("no molid matches with the HMDB identifier")

#TODO: hier kun je WHERE molid verwijderen.

def fetch_fragments_and_annotations_for_molid(molid: int,path_to_results_db_file: str) -> list:
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
    print(f"list with frag and smiles {list_with_fragments_and_smiles}")
    return list_with_fragments_and_smiles

def parse_line_with_motifs_and_querries(line: str) -> tuple:
    """
    Takes a line in a tab separated file and separates it into motif, feature list and massql query
    :param line: str, line in a tab separated file containing the motif, list of features and the massql query
    :return: returns the motif, feature_list and query in a tuple
    """
    motif=re.search(r'(.*)    (.*)    (.*)', line).group(1)
    feature_list=re.search(r'(.*)    (.*)    (.*)', line).group(2)
    query=re.search(r'(.*)    (.*)    (.*)', line).group(3)
    return motif,feature_list,query

def look_for_features(txt_file_with_motif_and_frag: str, current_motif: str) -> list:
    """
    Returns the list of features that make up the motif that are according to MassQL present in the spectrum file.

    :param txt_file_with_motif_and_frag: str, path to the output file from make_pdf_with_smiles.py in which each feature
    of each mass2motif is indicated.
    :param current_motif: str, motif in the name of the spectrum file, so the motif, which is according to MassQL,
    present in the spectrum.
    :return: list of feature belonging to the mass2motif which is according to MassQL present in the spectrum.
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
                # add the features of the current motif to the list of features 1 by 1.
                list_of_features.append(object[0])
            print(f"list features{list_of_features}")
            return list_of_features
    assert False, "The file with all motifs did not contain the specified motif"


def make_list_of_losses(list_with_fragments_and_smiles: list) -> list:
    """
    Makes a list of the Da of the losses between the parent mass and the annotated fragments

    :param list_with_fragments_and_smiles: list of tuples containing the fragment m/z and smiles for each fragment of
    the matched compound
    :return: a list with parent mass at position zero and after that the Da (Decimal) of the losses between the parent
    mass and the annotated fragments
    """
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
    print(f"list with losses {list_with_losses}")
    return list_with_losses

def get_smiles_of_loss(parent_string: str,fragment_string: str):
    """
    Takes the smiles of a parent ion and fragment ion and returns the smiles of the neutral loss

    :param parent_string: str, smiles of the precursor ion
    :param fragment_string: str, smiles of the annotated fragment
    :return: str, smiles of the neutral loss
    """
    parent = Chem.MolFromSmiles(parent_string)
    print(Chem.MolToSmiles(parent))
    fragment = Chem.MolFromSmiles(fragment_string)
    #print(Chem.MolToSmiles(fragment))
    # sometimes a parent or fragment string given by MAGMa cannot be converted into a real molecule and then an error
    # is thrown in this function, an ArgumentError... but how do I handle an error of this imported function?
    try:
        neutral_loss = Chem.ReplaceCore(parent, fragment)
    except:
        return None
    print(Chem.MolToSmiles(neutral_loss))
    # if the fragment smiles is not present in the parent smiles then neutral loss could be None
    if neutral_loss != None:
        try:
            neutral_loss = Chem.GetMolFrags(neutral_loss, asMols=True)
            print(neutral_loss)
        # the results from replace core could be a molecule that is not a real molecule upon which Chem.GetMolFrag will
        # give an error, so that function should be tried.
        except:
            print("neutral loss probably not a good molecule, so not added to results")
            return None
        print("this is neutral loss")
        print(neutral_loss)
        smiles_neutral_loss = ""
        for loss in neutral_loss:
            # if the neutral loss is a list of more substructures than the neutral loss is not 1 substructure, but a
            # splintered substructure and its hard to put back to a molecule that makes sense, so then its not returned.
            print(Chem.MolToSmiles(loss))
            part_of_smiles_loss = re.search(r'(\[.*\])(.*)',
                                        Chem.MolToSmiles(loss)).group(2)
            print(part_of_smiles_loss)
            # the neutral loss is added into one string
            smiles_neutral_loss+=part_of_smiles_loss
            if len(neutral_loss)>1:
                print("splintered substructure")
                try:
                    MolWt(Chem.MolFromSmiles(f'{smiles_neutral_loss}'))
                except:
                    return None
            print(smiles_neutral_loss)
        return smiles_neutral_loss

def search_for_smiles(list_of_features: list,list_with_fragments_and_smiles: list) -> list:
    """
    Makes a list of lists of features of the mass2motif that are annotated by MAGMa

    :param list_of_features: list of feature belonging to the mass2motif which is according to MassQL present in the
    spectrum.
    :param list_with_fragments_and_smiles: list of tuples containing the fragment m/z and smiles for each fragment of
    the matched compound
    :return: list_with_annotated_features: list of lists containing with each list containing the name of the feature,
    the smiles annotation by magma and the count which is 1.
    """
    # for every feature in the motif
    print(list_of_features)
    list_with_annotated_features = []
    for feature in list_of_features:
        # more than 1 feature could be annotated by MAGMa so make a list of lists with
        # if the feature in the motif is a loss
        if re.search(r'loss', feature) != None:
            # make a list of losses from the list_with_fragments_and_smiles from MAGMa
            list_of_losses=make_list_of_losses(list_with_fragments_and_smiles)
            # and for every loss in this list of losses, try to find a match to the loss of the feature
            for index,rounded_loss in enumerate(list_of_losses):
                # match the feature to the fragment/neutral loss on 2 decimals, rounded up or down.
                lower_bound=float(re.search(r'\_(.*\..{2}).*', feature).group(1))-0.01
                upper_bound=round(float(re.search(r'\_(.*)', feature).group(1)), 2)+0.01
                if float(rounded_loss) in np.arange(lower_bound, upper_bound+0.01, 0.01):
                    # then retrieve the corresponding fragment string and the parent string belonging to the loss
                    precusor_mz, precusor_smiles = list_with_fragments_and_smiles[0]
                    fragment_mz, fragment_smiles= list_with_fragments_and_smiles[index]
                    smiles_neutral_loss=get_smiles_of_loss(precusor_smiles, fragment_smiles)
                    if smiles_neutral_loss != None:
                        # check if the smiles has the same molecular mass as the loss reported of the feature
                        print(smiles_neutral_loss)
                        mol_weight_from_smiles=MolWt(Chem.MolFromSmiles(f'{smiles_neutral_loss}'))
                        print(mol_weight_from_smiles)
                        rounded_mol_weight_from_smiles=Decimal(mol_weight_from_smiles).quantize(Decimal('.1'),
                                                                          rounding=ROUND_DOWN)
                        print(rounded_mol_weight_from_smiles)
                        if rounded_mol_weight_from_smiles==float(re.search(r'\_(.*\..{1}).*', feature).group(1)):
                            count = 1
                            list_with_annotated_features.append([feature,smiles_neutral_loss, count])
                        else:
                            count = 1
                            list_with_annotated_features.append([feature,smiles_neutral_loss, count, mol_weight_from_smiles])
        # if feature is not a loss
        else:
            for fragment in list_with_fragments_and_smiles:
                print(fragment)
                fragment_mz, fragment_smiles = fragment
                rounded_fragment = Decimal(fragment_mz).quantize(Decimal('.01'),
                                                                          rounding=ROUND_DOWN)
                lower_bound = float(re.search(r'\_(.*\..{2}).*', feature).group(1)) - 0.01
                upper_bound = round(float(re.search(r'\_(.*)', feature).group(1)), 2) + 0.01
                print(np.arange(lower_bound, upper_bound+0.01, 0.01))
                print(float(rounded_fragment))
                if float(rounded_fragment) in np.arange(lower_bound, upper_bound+0.01, 0.01):
                    count=1
                    list_with_annotated_features.append([feature,fragment_smiles, count])
    print(list_with_annotated_features)
    return list_with_annotated_features

def make_output_file(path_to_txt_file_with_motif_and_frag: str) -> tuple:
    """
    Makes a copy of a tab-separated file, puts the file in a dataframe and removes the last column.

    :param path_to_txt_file_with_motif_and_frag: str, path to the output file from make_pdf_with_smiles.py which is a
    tab delimited text file with the selected motif, its fragments+probabilities and the MassQL query
    :return: tuple with the filepath of the output file and a dataframe with motifs and features
    """
    # This function should be executed one time and then you just continue adding stuff to dataframe
    path,filename = os.path.split(path_to_txt_file_with_motif_and_frag)
    file_path = Path(rf"{path}/motif_features_annotated.txt")
    if os.path.exists(file_path):
        assert False, f"The output file: {file_path} exists, remove it"
    shutil.copyfile(path_to_txt_file_with_motif_and_frag, file_path)
    df_with_motifs=pd.read_csv(file_path, sep="    ", engine='python', header=None)
    df_with_motifs.columns=["motif_name", "features", "massql_query"]
    df_with_motifs=df_with_motifs.set_index(["motif_name"])
    del df_with_motifs["massql_query"]
    # add a column where a list of lists of feature, annotation and count can go
    df_with_motifs["LoL_feature_annotation_counts"] = ""
    return os.path.abspath(file_path), df_with_motifs

def write_spectrum_output_to_df(list_with_annotated_features: list, df_with_motifs: pd.DataFrame,
                                    current_motif: str) -> pd.DataFrame:
    """
    Writes the obtained list of annotations for each feature to a dataframe.

    :param list_with_annotated_features: list of lists containing with each list containing the name of the feature,
    the smiles annotation by magma and the count which is 1.
    :param df_with_motifs: Pandas Dataframe with the motifs as an index and in one column a list
    of lists containing the feature, probability of the feature, the ratio of associated document with the feature,
    a list of the documents that contain the feature for every selected feature. The second column is a column to put
    the smiles annotations for the features
    :param current_motif: str, motif in the name of the spectrum file, so the motif, which is according to MassQL,
    present in the spectrum.
    :return: df_with_motifs, but the list_with_annotated_features have been added to the second column or the counts was
    increased if the same feature and annotated was already present in the lists of lists
    """
    for feature in list_with_annotated_features:
        cell=df_with_motifs.at[current_motif,"LoL_feature_annotation_counts"]
        # if there are already annotations in the cell see if they are similar to the current annotation
        if cell!="":
            for list_with_feature in cell:
                if list_with_feature[0] == feature[0]:
                    if list_with_feature[1] == feature[1]:
                        list_with_feature[2] = int(list_with_feature[2]) + 1
                    else:
                        cell.append(feature)
                else:
                    cell.append(feature)
        # if there are no annotations in the cell just add the list containing the feature, annotation and count.
        else:
            df_with_motifs.at[current_motif,"LoL_feature_annotation_counts"] = [(feature),]
    return df_with_motifs

def write_output_to_file(updated_df_with_motifs: pd.DataFrame, file_path: str) -> str:
    """Outputs a tab delimited file with the selected motifs, their features and the smiles annotations

    :param updated_df_with_motifs: Pandas Dataframe with the motifs as an index and in one column a list
    of lists containing the feature, probability of the feature, the ratio of associated document with the feature,
    a list of the documents that contain the feature for every selected feature. The second column contains a list of
    lists with each list containing the name of the feature, the SMILES annotation by MAGMa and the count.
    :param file_path: str, the filepath of the output file with the annotations which is in the same directory as
    the text file from make_pdf_with_smiles.py.
    :return: filepath where the results tab delimited file with the selected motifs, their features and the smiles
    annotations is stored.
    """
    # when you are at the end of the spectrum list you add the things to the dataframe
    output_file = open(file_path, "w")
    for index, row in updated_df_with_motifs.iterrows():
        output_file.write("{0}    {1}    {2}".format(index, updated_df_with_motifs.at[index, "features"],
                                                               updated_df_with_motifs.at[index, "LoL_feature_annotation_counts"]))
        output_file.write("\n")
    output_file.close()
    return os.path.abspath(file_path)

def main():
    #main function of the script
    # this script is made to return the smiles corresponding to the features of a motif of one spectrum_file_from_massql
    #step 0: parse input
    before_script=time.perf_counter()
    path_to_structures_database=argv[1]
    path_to_spectrum_files=argv[2]
    path_to_store_results_db=argv[3]
    path_to_txt_file_with_motif_and_frag=argv[4]
    path_to_spectrum_files = open(path_to_spectrum_files)
    # step 1: make a dataframe to put the annotations in for the features (so the output)
    file_path, df_with_motifs = make_output_file(path_to_txt_file_with_motif_and_frag)
    for line in path_to_spectrum_files:
        line = line.strip()
        path_to_spectrum_file = line.replace('\n', '')
        #step 2: parse name of spectrum file and initialize database to save results from MAGMa
        new_or_exists, path_to_results_db_file, identifier, current_motif=construct_path_to_db(path_to_spectrum_file, path_to_store_results_db)
        if new_or_exists=="new":
            initialize_db_to_save_results(path_to_results_db_file)
            print(path_to_results_db_file)
            after_init=time.perf_counter()
            print("init database {0}".format(after_init-before_script))
            # step 3: add spectrum to be annotated into the results database
            add_spectrum_into_db(path_to_results_db_file, path_to_spectrum_file, abs_intensity_thres=1000,
                                                 mz_precision_ppm=80, mz_precision_abs=0.01, spectrum_file_type='mgf',
                                                 ionisation=1)
            after_add_spectrum=time.perf_counter()
            print("add spectrum {0}".format(after_add_spectrum-after_init))
            #step 4: annotate spectrum with MAGMa and store output in results database
            annotate_spectrum_with_MAGMa(path_to_structures_database, path_to_results_db_file, max_num_break_bonds=10, structure_db="hmdb",
                              ncpus=1)
            after_annotate=time.perf_counter()
            print("annotate spectrum {0}".format(after_annotate-after_add_spectrum))
        # step 5: get the molid that has a HMDB identifier that corresponds to the HMDB identifier of the spectrum
        molid=get_molid_of_matched_compound(path_to_results_db_file, identifier)
        after_annotate=time.perf_counter()
        # It could be that MAGMa didn't find the correct compound then molid is going to be None
        if molid is not None:
            # step 6: get a list of features of the current motif
            list_of_features = look_for_features(path_to_txt_file_with_motif_and_frag, current_motif)
            after_look_for_motif = time.perf_counter()
            print("look for features {0}".format(after_look_for_motif - after_annotate))
            # step 7: get a list of the annotated fragments and smiles of the matched compound
            list_with_fragments_and_smiles=fetch_fragments_and_annotations_for_molid(molid, path_to_results_db_file)
            after_get_fragments = time.perf_counter()
            print("fetch fragments {0}".format(after_get_fragments - after_look_for_motif))
            # step 8: look for matches in the two lists so you end up with a list with only smiles for the features of the
            # current motif
            list_with_annotated_features=search_for_smiles(list_of_features,list_with_fragments_and_smiles)
            # The list with annotated features could be None if the features of the motif are not in the list of
            # annotated features by MAGMa or if for example the neutral loss is splintered across the molecule, and
            # not 1 side group.
            if list_with_annotated_features is not None:
                after_get_smiles = time.perf_counter()
                print("search smiles {0}".format(after_get_smiles - after_get_fragments))
                # step 9: adds the new annotations for the features from the spectrum file in the database
                df_with_motifs=write_spectrum_output_to_df(list_with_annotated_features, df_with_motifs, current_motif)
                print(df_with_motifs)
            else:
                print("none of the features could be annotated")
            print("one spectrum done")
    # step 10: writes the dataframe to a tab-delimited file
    write_output_to_file(df_with_motifs, file_path)

if __name__ == "__main__":
    main()