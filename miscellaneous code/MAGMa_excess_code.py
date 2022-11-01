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
    # sometimes a parent or fragment string given by MAGMa cannot be converted into a real molecule and then an error
    # is thrown in this function, an ArgumentError
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

def make_regex_for_loss(parent_string, fragment_string):
    new_fragment_string = re.escape(fragment_string)
    list_fragment_string = list(new_fragment_string)
    for index, letter in enumerate(list_fragment_string):
        if index == 0:
            if letter == "\\":
                letter_with_slash = "[\(|\)]*" + letter
                list_fragment_string[index] = letter_with_slash
            else:
                letter_with_slash = "[\(|\)]*" + letter + "[\(|\)]*"
                list_fragment_string[index] = letter_with_slash
        elif letter != "\\":
            if index != 0:
                letter_with_slash = letter + "[\(|\)]*"
                list_fragment_string[index] = letter_with_slash
    re_search_string = "".join(list_fragment_string)
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

def make_regex_for_loss_version2(parent_string, fragment_string):
        new_fragment_string = re.escape(fragment_string)
        list_fragment_string = list(new_fragment_string)
        for index, letter in enumerate(list_fragment_string):
            if index == 0:
                if letter == "\\":
                    letter_with_slash = "[\(|\)]*" + letter
                    list_fragment_string[index] = letter_with_slash
                else:
                    letter_with_slash = "[\(|\)]*" + letter + "[\(|\)]*"
                    list_fragment_string[index] = letter_with_slash
            elif letter != "\\":
                if index != 0:
                    letter_with_slash = letter + "[\(|\)]*"
                    list_fragment_string[index] = letter_with_slash
        re_search_string = "".join(list_fragment_string)
        if re.search(rf'(.+)({re_search_string})',
                     parent_string) != None and re.search(rf'({re_search_string})(.+)',
                                                          parent_string) != None:
            smiles_neutral_loss = re.search(rf'(.+)({re_search_string})',
                                            parent_string).group(1) + re.search(rf'({re_search_string})(.+)',
                                                                                parent_string).group(2)
            mols = [Chem.MolFromSmiles(parent_string), Chem.MolFromSmiles(smiles_neutral_loss)]
            print(mols)
            MCS_in_smart_string = rdFMCS.FindMCS(mols, ringMatchesRingOnly=True).smartsString
            print(MCS_in_smart_string)
            print(Chem.MolToSmiles(Chem.MolFromSmarts(MCS_in_smart_string)))
            return Chem.MolToSmiles(Chem.MolFromSmarts(MCS_in_smart_string))
        elif re.search(rf'(.+)({re_search_string})',
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