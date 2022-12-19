#!/usr/bin/env python3
"""
Author: Anna-Marie Seelen
Studentnumber:1008970
Description: This script runs MAGMa on one mgf formatted spectrum file and looks for smiles annotations for the m2m_frag_and_or_loss
of the corresponding mass2motif for which the spectrum file was made using MassQL
Usage: python3 MAGMa.py *path_to_mgf_spectrum_files_from_massql_script* *path_to_store_results_db*
*path_to_txt_file_with_massql_queries* *path_to_save_png_of_substruc*

    path_to_mgf_spectrum_files_from_massql_script: path where for each selected motif a file is stored with the
    mgf-style selected library spectra found by the MassQL script
    path_to_store_results_db: path where the output sqlite database from MAGMa with the annotations for the spectrum
    should be stored
    path_to_txt_file_with_massql_queries: path to the output file from make_pdf_with_smiles.py in which each fragment
    or neutral loss of each mass2motif is indicated.
    path_to_save_png_of_substruc: path to where the structures of all the library spectra with the annotated Mass2Motif
    fragment highlighted can be stored
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
from rdkit.Chem.Descriptors import MolWt
import time
import shutil
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import DrawingOptions, MolDrawOptions
import signal

# functions
def parse_input(mgf_file: str) -> dict:
    """Parses mgf into strings, where each string is a spectrum and stores those in dict with the spectrum_id as key.

    :param mgf_file: str, name of mgf formatted file containing MS/MS spectra from GNPS
    :return: dictionary with {spectrum_id:record} where each record is a string containing the mgf-style accession of
    one compound
    """

    lines_mgf_file = open(mgf_file)
    spectrum_record_bool = False
    mgf_spectrum_records = []
    mgf_spectrum_record = ""
    for line in lines_mgf_file:
        if line.startswith("BEGIN IONS"):
            spectrum_record_bool = True
        if spectrum_record_bool:
            mgf_spectrum_record += line
        if line.startswith("END IONS"):
            mgf_spectrum_records.append(mgf_spectrum_record)
            spectrum_record_bool = False
            mgf_spectrum_record = ""

    dict_with_mgf_spectra = {}
    for mgf_spectrum_record in mgf_spectrum_records:
        mgf_spectrum_record = mgf_spectrum_record.strip()
        # look for the spectrumid in the string and use it as a key for the dict
        key = re.search(r'SPECTRUMID=(.*)', mgf_spectrum_record).group(1)
        if key is not None:
            dict_with_mgf_spectra[key] = mgf_spectrum_record
    return dict_with_mgf_spectra

def construct_path_to_db(spectrum_id, path_to_store_results_db: str) -> tuple:
    """

    :param spectrum_id: str, the spectrum_id of GNPS spectra files formatted like this: CCMSLIB00000425029
    :param path_to_store_results_db: str, path where the output sqlite database from MAGMa with the annotations for the
    library spectrum should be stored
    :return: tuple, first element is a string which indicates if the database already exists or not, second element
    is the path of the results database file of one library spectrum
    """
    # generate a good database name with the spectrum id in it.
    file_path_results_db = Path(fr"{path_to_store_results_db}/MAGMa_db_{spectrum_id}.sqlite")
    # if you try to annotate something in an existing database it will go wrong, so if the database exists do not
    # annotate it again. Here you assume that all the magma steps init, add_spectrum, add_smiles, annotate have been
    # done
    if os.path.exists(file_path_results_db):
        # assert False, "The results db for this spectrum already exists, remove it"
        return "exists", file_path_results_db
    else:
        return "new", file_path_results_db

def make_mgf_txt_file_for_spectrum(spectrum_id: str, spectrum_record_mgf: str, path_to_store_spectrum_files: str) -> str:
    """
    Writes one spectrum containing the motif to a text file in mgf format using the identifier given by MassQL

    :param spectrum_id: str, the spectrum_id of GNPS spectra files formatted like this: CCMSLIB00000425029
    :param spectrum_record_mgf: str, a string with the mgf-formatted spectrum of the current spectrum_id
    :param path_to_store_spectrum_files: str, folder where all the mgf-formatted text files with spectra will be stored.
    :return: path of file with one mgf-formatted library spectrum with the GNPS identifier of the spectrum in the file
    name.
    """

    # create the name for the spectrum file with the GNPS identifier
    path_to_spectrum_file = Path(fr"{path_to_store_spectrum_files}/spectrum_{spectrum_id}.txt")
    # write the mgf spectrum saved in the dict only to a text file if the text file doesn't exist already, if the text
    # file does exist it is assumed to be the right text file in the right formatting
    if os.path.exists(path_to_spectrum_file):
        return os.path.abspath(path_to_spectrum_file)
    else:
        spectrum_file = open(path_to_spectrum_file, "w")
        spectrum_file.write(spectrum_record_mgf)
        spectrum_file.close()
        return os.path.abspath(path_to_spectrum_file)

def initialize_db_to_save_results(path_to_results_db: str) -> None:
    """
    Initializes a sqlite database using the MAGMa CL init_db function to store annotation results from MAGMa for spectrum

    :param path_to_store_results_db: str, path where the output database from MAGMa with the annotations for the
    spectrum should be stored.
    :return: None
    """

    cmd = f'magma init_db {path_to_results_db}'
    try:
        e = subprocess.check_call(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    # For some reason the text that MAGMa returns on the command line is seen as an error so include this line to keep
    # script from stopping.
    except subprocess.CalledProcessError:
        pass
    return None

def add_spectrum_into_db(path_to_results_db_file: str, path_to_spectrum_file: str,
                         spectrum_file_type='mgf', ionisation=1) -> None:
    """
    Adds the spectrum to be annotated in the sqlite results database using the MAGMa CL read_ms_data function.

    :param path_to_results_db_file: str, path to the results database in which the library spectrum and peak annotations
    will be stored
    :param path_to_spectrum_file: str, path where the spectrum file is stored with a spectrum from the HMDB that
    contains a motif which is determined with MassQL.py script.
    :param spectrum_file_type: str, indicates how the spectrum_file is formatted (default: mgf)
    :param ionisation: {1, -1}, indicates the ionisation mode, -1 is negative, 1 is positive (default: 1)
    :return: None
    """
    cmd = f'magma read_ms_data -f {spectrum_file_type} -i {ionisation} {path_to_spectrum_file} {path_to_results_db_file}'
    # sometimes read_ms_data gives an error that is not usefull at all so:
    try:
        e = subprocess.check_call(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError:
        print(f"read_ms data gave an error, so this sqlite lib was not made for {path_to_spectrum_file}")
        os.remove(path_to_results_db_file)
        pass
    return None

def add_structure_to_db(path_to_results_db_file: str, smiles: str) -> None:
    """
    Adds the smiles of the precursor molecule to the results database

    :param path_to_results_db_file: str, path to the results database in which the library spectrum and peak annotations
    will be stored
    :param smiles: the smiles of the precursor molecule belonging to the library MS/MS spectrum
    :return: None
    """
    cmd = f"""magma add_structures -t smiles '{smiles}' {path_to_results_db_file}"""
    e = subprocess.check_call(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    return None

def annotate_spectrum_with_MAGMa(path_to_results_db_file: str,
                                 max_num_break_bonds=10, ncpus=1) -> None:
    """
    Annotates generated fragments using MAGMa CL annotate function, which uses the precursor smiles.

    :param path_to_results_db_file: str, path to the results database in which the library spectrum and peak annotations
    will be stored
    :param max_num_break_bonds: int, maximum number of bond breaks to generate substructures (default: 10)
    :param ncpus: int, number of parallel cpus to use for annotation (default: 1)
    :return: None
    """
    cmd = f"magma annotate -b {max_num_break_bonds} -n {ncpus} {path_to_results_db_file}"
    e = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
    return None

def fetch_fragments_and_annotations(path_to_results_db_file: str) -> list:
    """
    Makes list of tuples containing the fragment m/z and smiles of each fragment peak

    :param path_to_results_db_file: str, path to the results database in which the library spectrum and peak annotations
    are stored
    :return: list of tuples containing the fragment m/z and annotated smiles for each fragment peak
    """
    conn = sqlite3.connect(path_to_results_db_file)
    sqlite_command = \
        f"""SELECT mz,smiles FROM fragments"""
    cur = conn.cursor()
    cur.execute(sqlite_command)
    list_with_fragments_and_smiles = cur.fetchall()
    return list_with_fragments_and_smiles

def parse_line_with_motif_and_query(line: str) -> tuple:
    """
    Takes a line in a tab separated file and separates it into motif, m2m_frag_or_loss list and massql query
    :param line: str, line in a tab separated file containing the motif, list of m2m_frag_and_or_loss and the massql query
    :return: returns the motif, m2m_frag_or_loss_list and the massql query in a tuple
    """
    splitted_line = line.split("\t")
    assert len(
        splitted_line) == 3, "Expected a file with lines with tabs separating motif, m2m_frag_or_loss_list and massql query"
    motif, m2m_frag_and_or_loss, massql_query = splitted_line
    return motif, m2m_frag_and_or_loss, massql_query

def rearrange_m2m_frag_and_or_loss(m2m_frag_and_or_loss_list: str) -> list:
    """
    Returns the list of m2m_frag_and_or_losses that make up the motif that are according to MassQL present in the spectrum file.

    :param m2m_frag_and_or_loss_list: list of the
    #TODO: what is this?
    :return: list of m2m_frag_or_loss belonging to the mass2motif which is according to MassQL present in the spectrum.
    """
    m2m_frag_and_or_loss_list_1 = ast.literal_eval(m2m_frag_and_or_loss_list)
    list_of_m2m_frag_and_or_loss = []
    for object in m2m_frag_and_or_loss_list_1:
        # add the m2m_frag_and_or_loss of the motif to the list of m2m_frag_and_or_loss 1 by 1.
        list_of_m2m_frag_and_or_loss.append(object[0])
    return list_of_m2m_frag_and_or_loss

def make_list_of_losses(list_with_fragments_and_smiles: list) -> list:
    """
    Makes a list of the Da of the losses between the parent mass and the annotated fragments

    :param list_with_fragments_and_smiles: list of tuples containing the fragment m/z and smiles for each fragment of
    the matched compound
    :return: a list with parent mass at position zero and after that the Da (Decimal) of the losses between the parent
    mass and the annotated fragments
    """
    list_with_losses = []
    for index, fragment in enumerate(list_with_fragments_and_smiles):
        if index != 0:
            precusor_mz, precusor_smiles = list_with_fragments_and_smiles[0]
            fragment_mz, fragment_smiles = fragment
            calculated_loss = Decimal(precusor_mz - fragment_mz).quantize(Decimal('.01'),
                                                                          rounding=ROUND_DOWN)
            list_with_losses.append(calculated_loss)
        else:
            precusor_mz, precusor_smiles = fragment
            list_with_losses.append(precusor_mz)
    return list_with_losses


def get_mol_block(path_to_results_db_file: str) -> str:
    """
    Retrieves the molblock of the precursor molecule from the results database from MAGMa

    :param path_to_results_db_file: str, path to the results database in which the spectrum and annotations are stored.
    :return: the molblock in string format of the precursor molecule representing in numbers what the mol looks like
    """
    conn = sqlite3.connect(path_to_results_db_file)
    sqlite_command = \
        f"""SELECT mol FROM molecules"""
    cur = conn.cursor()
    cur.execute(sqlite_command)
    molblock = cur.fetchall()
    # The mol block is in a list of tuples
    molblock = molblock[0][0]
    return molblock

def get_atom_list(path_to_results_db_file: str, mass_of_frag: float) -> list:
    """
    Retrieves the atom_list for a fragment peak in the results database from magma with a particular mass

    :param path_to_results_db_file: str, path to the results database in which the spectrum and annotations are stored.
    :param mass_of_frag: float, the mass of the fragment for which a atom_list should be retrieved
    :return: the atom list which is a list of integers representing the atoms of the precursor mol that the
    fragment contains
    """
    conn = sqlite3.connect(path_to_results_db_file)
    sqlite_command = \
        f"""SELECT atoms FROM fragments WHERE mz = {mass_of_frag}"""
    cur = conn.cursor()
    cur.execute(sqlite_command)
    atom_list = cur.fetchall()
    # The atom list is in a list of tuples
    atom_list = atom_list[0][0]
    atoms = [int(a) for a in atom_list.split(',')]
    return atoms

def get_bond_list(atom_list, molblock):
    """

    :param atom_list: list, the atom list which is a list of integers representing the atoms of the precursor mol that
    the fragment contains
    :param molblock: str, molblock in string format of the precursor molecule representing in numbers what the mol looks
    like
    :return: list, list of int representing the bonds between the atoms of the fragment
    """
    # creating a bond list for visualization
    mol = Chem.MolFromMolBlock(molblock)

    # creating a bond list for visualization
    bond_list = []
    bonds_in_prec_mol = [(x.GetBeginAtomIdx(), x.GetEndAtomIdx()) for x in mol.GetBonds()]
    for bond in bonds_in_prec_mol:
        atom_1,atom_2=bond
        #if atom 1 and atom 2 are also both in the atom_list then they are connected to each other in the molecule!
        if all(x in atom_list for x in [atom_1,atom_2]):
            bond_list.append(mol.GetBondBetweenAtoms(atom_1, atom_2).GetIdx())
    return bond_list

def loss2smiles(molblock, atomlist):
    """
    Get the atom list and bond list and the rdkit.Mol class of the neutral loss

    :param molblock: str, molblock in string format of the precursor molecule representing in numbers what the mol looks
    like
    :param atomlist: list, the atom list which is a list of integers representing the atoms of the precursor mol that
    the fragment contains. This is the atom_list of the fragment that when subtracted from the precursor mol is the
    neutral loss we are looking for
    :return tuple with the atom_list, the bond_list and the rdkit.Mol of the neutral loss
    """
    # getting the neutral loss in mol object form
    mol = Chem.MolFromMolBlock(molblock)
    emol = Chem.EditableMol(mol)
    neutral_loss_atom_list=list(reversed(range(mol.GetNumAtoms())))

    for atom in reversed(range(mol.GetNumAtoms())):
        if atom in atomlist:
            neutral_loss_atom_list.remove(atom)
            emol.RemoveAtom(atom)
    neutral_loss = emol.GetMol()

    # creating a bond list for visualization
    neutral_loss_bond_list = []
    bonds_in_prec_mol = [(x.GetBeginAtomIdx(), x.GetEndAtomIdx()) for x in mol.GetBonds()]
    for bond in bonds_in_prec_mol:
        atom_1,atom_2=bond
        #if atom 1 and atom 2 are also both in the atom_list then they are connected to each other in the molecule!
        if all(x in neutral_loss_atom_list for x in [atom_1,atom_2]):
            neutral_loss_bond_list.append(mol.GetBondBetweenAtoms(atom_1, atom_2).GetIdx())

    return neutral_loss_atom_list, neutral_loss_bond_list, Chem.MolToSmiles(neutral_loss)


def search_for_smiles(list_of_m2m_frag_and_or_loss: list, list_with_fragments_and_smiles: list, path_to_results_db_file: str, spectrum_id: str,
                      motif: str, path_to_save_png_of_substruc: str) -> list:
    """
    Makes a list of lists of m2m_frag_and_or_loss of the mass2motif that are annotated by MAGMa

    :param list_of_m2m_frag_and_or_loss: list of fragment and/or loss belonging to the mass2motif which is according to
    MassQL present in the library spectrum.
    :param list_with_fragments_and_smiles: list of tuples containing the fragment m/z and smiles for each annotated
    fragment peak
    :param path_to_results_db_file: str, path to the results database in which the library spectrum and peak annotations
    are stored
    :param spectrum_id: str, the spectrum_id of GNPS spectra files formatted like this: CCMSLIB00000425029
    :param motif: str, name of the current motif
    :param path_to_save_png_of_substruc: path to where the structures of all the library spectra with the annotated
    Mass2Motif fragment highlighted can be stored
    :return: list_with_annotated_m2m_frag_and_or_loss: list of lists containing with each list containing the name of
    the m2m_frag_or_loss, the smiles annotation by magma and the count which is 1.
    """
    # for every m2m_frag_or_loss in the motif
    list_with_annotated_m2m_frag_and_or_loss = []
    for m2m_frag_or_loss in list_of_m2m_frag_and_or_loss:
        # more than 1 m2m_frag_or_loss could be annotated by MAGMa so make a list of lists with
        # if the m2m_frag_or_loss in the motif is a loss
        if re.search(r'loss', m2m_frag_or_loss) != None:
            # make a list of losses from the list_with_fragments_and_smiles from MAGMa
            list_of_losses = make_list_of_losses(list_with_fragments_and_smiles)
            # and for every loss in this list of losses, try to find a match to the loss of the m2m_frag_or_loss
            for index, rounded_loss in enumerate(list_of_losses):
                # match the m2m_frag_or_loss to the fragment/neutral loss on 2 decimals, rounded up or down.
                lower_bound = float(re.search(r'\_(.*\..{2}).*', m2m_frag_or_loss).group(1)) - 0.01
                upper_bound = round(float(re.search(r'\_(.*)', m2m_frag_or_loss).group(1)), 2) + 0.01
                if float(rounded_loss) in np.arange(lower_bound, upper_bound + 0.01, 0.01):
                    # then retrieve the corresponding fragment string and the parent string belonging to the loss
                    precursor_mz, precursor_smiles = list_with_fragments_and_smiles[0]
                    fragment_mz, fragment_smiles = list_with_fragments_and_smiles[index]
                    molblock = get_mol_block(path_to_results_db_file)
                    atomlist = get_atom_list(path_to_results_db_file, float(fragment_mz))
                    neutral_loss_atom_list, neutral_loss_bond_list, smiles_neutral_loss = loss2smiles(molblock,
                                                                                                      atomlist)
                    # if the smiles of the neutral_loss could be found with the loss2smiles function
                    if smiles_neutral_loss != None:
                        # check if the smiles has the same molecular mass as the loss reported of the m2m_frag_or_loss
                        vis_substructure_in_precursor_mol(path_to_save_png_of_substruc, molblock, neutral_loss_atom_list,
                                                          neutral_loss_bond_list, spectrum_id, motif)
                        # sometimes the molecular weight cannot be calculated if the loss is from a cyclic molecule
                        # because some atoms will be lowercase, but the molecule will not be aromatic.
                        try:
                            mol_weight_from_smiles = MolWt(Chem.MolFromSmiles(f'{smiles_neutral_loss}'))
                        # if the molecular weight cannot be calculated convert all the lowercase letters to uppercase,
                        # because then you will get a molecule that is not aromatic and does not have aromatic atoms
                        except:
                            all_uppercase_smiles = smiles_neutral_loss.upper()
                            mol_weight_from_smiles = MolWt(Chem.MolFromSmiles(f'{all_uppercase_smiles}'))

                        rounded_mol_weight_from_smiles = Decimal(mol_weight_from_smiles).quantize(Decimal('.1'),
                                                                                                  rounding=ROUND_DOWN)
                        if rounded_mol_weight_from_smiles == float(re.search(r'\_(.*\..{1}).*', m2m_frag_or_loss).group(1)):
                            count = 1
                            list_with_annotated_m2m_frag_and_or_loss.append([m2m_frag_or_loss, smiles_neutral_loss, count])
                        # if the mass is not the same as the reported m2m_frag_or_loss then the mass is reported in the output file
                        else:
                            count = 1
                            list_with_annotated_m2m_frag_and_or_loss.append(
                                [m2m_frag_or_loss, smiles_neutral_loss, count, mol_weight_from_smiles])

        # if m2m_frag_or_loss is not a loss
        else:
            for fragment in list_with_fragments_and_smiles:
                fragment_mz, fragment_smiles = fragment
                rounded_fragment = Decimal(fragment_mz).quantize(Decimal('.01'),
                                                                 rounding=ROUND_DOWN)
                lower_bound = float(re.search(r'\_(.*\..{2}).*', m2m_frag_or_loss).group(1)) - 0.01
                upper_bound = round(float(re.search(r'\_(.*)', m2m_frag_or_loss).group(1)), 2) + 0.01
                if float(rounded_fragment) in np.arange(lower_bound, upper_bound + 0.01, 0.01):
                    atomlist = get_atom_list(path_to_results_db_file, float(fragment_mz))
                    precursor_mz, precursor_smiles = list_with_fragments_and_smiles[0]
                    if fragment_smiles != None:
                        molblock = get_mol_block(path_to_results_db_file)
                        bond_list = get_bond_list(atomlist, molblock)
                        vis_substructure_in_precursor_mol(path_to_save_png_of_substruc, molblock, atomlist, bond_list, spectrum_id, motif)
                    count = 1
                    list_with_annotated_m2m_frag_and_or_loss.append([m2m_frag_or_loss, fragment_smiles, count])
    return list_with_annotated_m2m_frag_and_or_loss


def vis_substructure_in_precursor_mol(path_to_save_png_of_substruc: str,mol_block: str, atom_list: list, bond_list: list, spectrum_id: str, motif: str) -> None:
    """
    Make a png image file for the library structure and highlighted the frag or loss annotation according to MAGMa

    :param path_to_save_png_of_substruc: path to where the structures of all the library spectra with the annotated Mass2Motif
    fragment highlighted can be stored
    :param mol_block: str, a molblock of the precursor molecule representing in numbers what the structure looks like
    :param atom_list: list, the atom list which is a list of integers representing the atoms of the precursor mol that
    the fragment or loss contains.
    :param bond_list: list, list of int representing the bonds between the atoms of the fragment
    :param spectrum_id: str, the spectrum_id of GNPS spectra files formatted like this: CCMSLIB00000425029
    :param motif: str, name of the Mass2Motif for which this annotation was obtained
    :return: None
    """
    opts = MolDrawOptions()
    # make the structure black, so no colours for oxygen atom and stuff
    opts.updateAtomPalette({k: (0, 0, 0) for k in DrawingOptions.elemDict.keys()})

    mol = Chem.MolFromMolBlock(mol_block)
    # the highlight colour is still red in the resulting png although it is different in this function.
    # So changing the highlighting colour in this function doesn't seem to work
    Draw.MolToFile(mol, f"{path_to_save_png_of_substruc}/{spectrum_id}_{motif}.png",
                   highlightAtoms=atom_list, highlightBonds=bond_list, highlightColor=[0,0,0], options=opts)
    return None

def make_output_file(path_to_txt_file_with_massql_queries: str) -> tuple:
    """
    Makes a copy of a tab-separated file, puts the file in a dataframe and removes the last column.

    :param path_to_txt_file_with_massql_queries: str, path to the output file from make_pdf_with_smiles.py which is a
    tab delimited text file with the selected motif, its fragments+probabilities and the MassQL query
    :return: tuple with the filepath of the output file and a dataframe with motifs and m2m_frag_and_or_loss
    """
    # This function should be executed one time and then you just continue adding stuff to dataframe
    path, filename = os.path.split(path_to_txt_file_with_massql_queries)
    file_path = Path(rf"{path}/motif_m2m_frag_and_or_loss_annotated.txt")
    if os.path.exists(file_path):
        assert False, f"The output file: {file_path} exists, remove it"
    shutil.copyfile(path_to_txt_file_with_massql_queries, file_path)
    df_with_motifs = pd.read_csv(file_path, sep="\t", engine='python', header=None)
    df_with_motifs.columns = ["motif_name", "m2m_frag_and_or_loss", "massql_query"]
    df_with_motifs = df_with_motifs.set_index(["motif_name"])
    del df_with_motifs["massql_query"]
    # add a column where a list of lists of m2m_frag_or_loss, annotation and count can go
    df_with_motifs["LoL_m2m_frag_or_loss_annotation_counts"] = ""
    return os.path.abspath(file_path), df_with_motifs

def write_spectrum_output_to_df(list_with_annotated_m2m_frag_and_or_loss: list, df_with_motifs: pd.DataFrame,
                                current_motif: str) -> pd.DataFrame:
    """
    Writes the obtained list of annotations for each m2m_frag_or_loss to a dataframe.

    :param list_with_annotated_m2m_frag_and_or_loss: list of lists containing with each list containing the name of the m2m_frag_or_loss,
    the smiles annotation by magma and the count which is 1.
    :param df_with_motifs: Pandas Dataframe with the motifs as an index and in one column a list
    of lists containing the m2m_frag_or_loss, probability of the m2m_frag_or_loss, the ratio of associated document with the m2m_frag_or_loss,
    a list of the documents that contain the m2m_frag_or_loss for every selected m2m_frag_or_loss. The second column is a column to put
    the smiles annotations for the m2m_frag_and_or_loss
    :param current_motif: str, motif in the name of the spectrum file, so the motif, which is according to MassQL,
    present in the spectrum.
    :return: df_with_motifs, but the list_with_annotated_m2m_frag_and_or_loss have been added to the second column or the counts was
    increased if the same m2m_frag_or_loss and annotated was already present in the lists of lists
    """
    for m2m_frag_or_loss in list_with_annotated_m2m_frag_and_or_loss:
        # if there are already annotations in the cell see if they are similar to the current annotation
        if df_with_motifs.at[current_motif, "LoL_m2m_frag_or_loss_annotation_counts"] != "":
            list_of_m2m_frag_and_or_loss_in_df = [list_with_m2m_frag_or_loss[0] for list_with_m2m_frag_or_loss in df_with_motifs.at[current_motif, "LoL_m2m_frag_or_loss_annotation_counts"]]
            if m2m_frag_or_loss[0] not in list_of_m2m_frag_and_or_loss_in_df:
                df_with_motifs.at[current_motif, "LoL_m2m_frag_or_loss_annotation_counts"].append(m2m_frag_or_loss)
            else:
                list_of_smiles_in_df = [list_with_m2m_frag_or_loss[1] for list_with_m2m_frag_or_loss in df_with_motifs.at[current_motif, "LoL_m2m_frag_or_loss_annotation_counts"]]
                if m2m_frag_or_loss[1] not in list_of_smiles_in_df:
                    df_with_motifs.at[current_motif, "LoL_m2m_frag_or_loss_annotation_counts"].append(m2m_frag_or_loss)
                else:
                    for list_with_m2m_frag_or_loss in df_with_motifs.at[current_motif, "LoL_m2m_frag_or_loss_annotation_counts"]:
                        if list_with_m2m_frag_or_loss[0] == m2m_frag_or_loss[0]:
                            if list_with_m2m_frag_or_loss[1] == m2m_frag_or_loss[1]:
                                list_with_m2m_frag_or_loss[2] = int(list_with_m2m_frag_or_loss[2]) + int(m2m_frag_or_loss[2])

        # if there are no annotations in the cell just add the list containing the m2m_frag_or_loss, annotation and count.
        else:
            df_with_motifs.at[current_motif, "LoL_m2m_frag_or_loss_annotation_counts"] = [m2m_frag_or_loss, ]
    df_with_motifs.at[current_motif, "LoL_m2m_frag_or_loss_annotation_counts"]=sorted(df_with_motifs.at[current_motif, "LoL_m2m_frag_or_loss_annotation_counts"], key=lambda x: x[2],
           reverse=True)
    return df_with_motifs

def write_output_to_file(updated_df_with_motifs: pd.DataFrame, file_path: str) -> str:
    """Outputs a tab delimited file with the selected motifs, their m2m_frag_and_or_loss and the smiles annotations

    :param updated_df_with_motifs: Pandas Dataframe with the motifs as an index and in one column a list
    of lists containing the m2m_frag_or_loss, probability of the m2m_frag_or_loss, the ratio of associated document with the m2m_frag_or_loss,
    a list of the documents that contain the m2m_frag_or_loss for every selected m2m_frag_or_loss. The second column contains a list of
    lists with each list containing the name of the m2m_frag_or_loss, the SMILES annotation by MAGMa and the count.
    :param file_path: str, the filepath of the output file with the annotations which is in the same directory as
    the text file from make_pdf_with_smiles.py.
    :return: filepath where the results tab delimited file with the selected motifs, their m2m_frag_and_or_loss and the smiles
    annotations is stored.
    """
    # when you are at the end of the spectrum list you add the things to the dataframe
    output_file = open(file_path, "w")
    for index, row in updated_df_with_motifs.iterrows():
        output_file.write("{0}    {1}    {2}".format(index, updated_df_with_motifs.at[index, "m2m_frag_and_or_loss"],
                                                     updated_df_with_motifs.at[index, "LoL_m2m_frag_or_loss_annotation_counts"]))
        output_file.write("\n")
    output_file.close()
    return os.path.abspath(file_path)

class timeout:
    """
    Class is copied from internet, because sometimes MAGMa annotate is really slow and doesn't annotate. This
    Class prevents magma annotate from running too long.
    """
    def __init__(self, seconds=1, error_message='Timeout'):
        self.seconds = seconds
        self.error_message = error_message

    def handle_timeout(self, signum, frame):
        raise TimeoutError(self.error_message)

    def __enter__(self):
        signal.signal(signal.SIGALRM, self.handle_timeout)
        signal.alarm(self.seconds)

    def __exit__(self, type, value, traceback):
        signal.alarm(0)

def main():
    # main function of the script
    # step 0: parse input
    before_script = time.perf_counter()
    path_to_store_spectrum_files = argv[1]
    path_to_store_results_db = argv[2]
    path_to_txt_file_with_massql_queries = argv[3]
    path_to_save_png_of_substruc = argv[4]
    # step 1: make a dataframe to put the annotations in for the m2m_frag_and_or_loss (so the output)
    file_path, df_with_motifs = make_output_file(path_to_txt_file_with_massql_queries)
    with open(path_to_txt_file_with_massql_queries, "r") as lines_motif_file:
        for line in lines_motif_file:
            line = line.strip()
            line = line.replace('\n', '')
            # step 2: retrieve the motif and the massql query one by one
            motif, m2m_frag_and_or_loss, massql_query = parse_line_with_motif_and_query(line)
            assert os.path.exists(
                fr'{path_to_store_spectrum_files}/mgf_spectra_for_{motif}_from_massql.txt'), f"there is no mgf file with spectra to annotate for {motif}"
            # for the spectra files from MassQL look for the file with the mgf spectra for the current motif
            for file in os.listdir(path_to_store_spectrum_files):
                if re.findall(fr'mgf_spectra_for_{motif}_from_massql.txt', file) != None:
                    if len(re.findall(fr'mgf_spectra_for_{motif}_from_massql.txt', file)) > 1:
                        assert False, f"more than 1 file with mgf style spectra for {motif} in dir {path_to_store_spectrum_files}"
            amount_of_annotated_spectra = 0
            # step 3: parse the mgf library spectra for the mass2motif into a dictionary
            dict_with_mgf_spectra = parse_input(f"{path_to_store_spectrum_files}/mgf_spectra_for_{motif}_from_massql.txt")
            for spectrum_id in dict_with_mgf_spectra.keys():
                # step 4: construct a results database path for a library spectrum for a mass2motif
                new_or_exists, path_to_results_db_file = construct_path_to_db(spectrum_id,
                                                                                      path_to_store_results_db)
                # if the database path with the spectrum identifier is not present then make the database from scratch
                if new_or_exists == "new":
                    # step 5: initialize the results database in which MAGMa will store the annotated fragment peaks
                    initialize_db_to_save_results(path_to_results_db_file)
                    # get the mgf style spectrum from the dictionary with mgf spectra
                    mgf_spectrum_record = dict_with_mgf_spectra[spectrum_id]
                    # step 6: write the current spectrum to a separate text file
                    path_to_spectrum_file = make_mgf_txt_file_for_spectrum(spectrum_id, mgf_spectrum_record,
                                                                                   path_to_store_spectrum_files)
                    # step 7: add spectrum to be annotated into the results database
                    add_spectrum_into_db(path_to_results_db_file, path_to_spectrum_file,
                                                         spectrum_file_type='mgf',
                                                         ionisation=1)
                    if os.path.exists(path_to_results_db_file):
                        # search the smiles in the spectrum record
                        smiles = re.search(r'SMILES=(.*)', mgf_spectrum_record).group(1)
                        # step 8: add the smiles of the spectrum to be annotated to the results database of massql
                        add_structure_to_db(path_to_results_db_file, smiles)
                        # remove the spectrum file with a single mgf-style library spectrum
                        os.remove(path_to_spectrum_file)
                    else:
                        return None
                    # step 9: annotate spectrum with MAGMa and store output in results database
                    # sometimes the annotate function takes WAY too long, so if it takes longer than 2 min it is skipped
                    with timeout(seconds=120):
                        try:
                            annotate_spectrum_with_MAGMa(path_to_results_db_file,
                                                         max_num_break_bonds=10,
                                                         ncpus=1)
                        except TimeoutError:
                            # os.remove(path_to_results_db_file)
                            continue
                # step 10: get a list of m2m_frag_and_or_loss of the current motif
                list_of_m2m_frag_and_or_loss = rearrange_m2m_frag_and_or_loss(m2m_frag_and_or_loss)
                # step 11: get a list of the annotated fragments and smiles of the annotated library spectrum
                list_with_fragments_and_smiles = fetch_fragments_and_annotations(path_to_results_db_file)
                # step 12: look for matches in the two lists so you end up with a list with only smiles for the m2m_frag_and_or_loss of the
                # current motif
                list_with_annotated_m2m_frag_and_or_loss = search_for_smiles(list_of_m2m_frag_and_or_loss,
                                                                         list_with_fragments_and_smiles,
                                                                         path_to_results_db_file, spectrum_id, motif, path_to_save_png_of_substruc)
                # The list with annotated m2m_frag_and_or_loss could be None if the m2m_frag_and_or_loss of the motif are not in the list of
                # annotated m2m_frag_and_or_loss by MAGMa: which means the Mass2Motif fragment or loss could not be annotated by MAGMa
                if list_with_annotated_m2m_frag_and_or_loss:
                    # step 13: add the new annotations for the m2m_frag_and_or_loss from the spectrum file in the database
                    amount_of_annotated_spectra += 1
                    df_with_motifs = write_spectrum_output_to_df(list_with_annotated_m2m_frag_and_or_loss, df_with_motifs,motif)
            print(f"the amount of spectra found in MassQL for motif {motif}: {len(dict_with_mgf_spectra.keys())}")
            print(f"amount of spectra that had an annotation for {motif} is {amount_of_annotated_spectra}")
    # step 14: writes the dataframe to a tab-delimited file
    write_output_to_file(df_with_motifs, file_path)
    after_script = time.perf_counter()
    print("how long the total script took {0}".format(after_script - before_script))

if __name__ == "__main__":
    main()