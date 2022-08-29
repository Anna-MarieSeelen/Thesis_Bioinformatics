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
from massql import msql_engine
from sys import argv
import sys
#TODO: fix that it imports the file, because the massql function sucksss
#import ALL_GNPS_210409_positive_processed_annotated_CF_NPC_classes.txt
import ntpath
import pyarrow.feather as feather
import os
import pyteomics
import re
# import rdkit.Chem as Chem
# from rdkit.Chem.Draw import MolToImage
#from fpdf import FPDF
import json
from pandas import json_normalize
import pandas as pd
import ast
import subprocess
from pathlib import Path
import re
import os, glob
import matchms
from matchms import Scores, Spectrum

# functions

def parse_line_with_motifs_and_querries(line):
    motif=re.search(r'(.*)    (.*)    (.*)', line).group(1)
    fragments=re.search(r'(.*)    (.*)    (.*)', line).group(2)
    query=re.search(r'(.*)    (.*)    (.*)', line).group(3)
    return motif,fragments,query

def make_json_file(pickle_file, path_to_store_json_file):
    obj = pd.read_pickle(pickle_file)
    spectrum_list=[]
    for spectrum in obj:
        spectrum_list.append(spectrum)
    matchms.exporting.save_as_json(spectrum_list,path_to_store_json_file)
    return None

def try_massql(query, json_file):
    df=msql_engine.process_query(query,json_file)
    df.rename(columns={'scan': 'spectrum_id'}, inplace=True)
    df.set_index("spectrum_id", inplace=True, drop=True)
    return df

def parse_input(filename):
    """Parses argonaut formatted file to extract accession number, organism name and DNA_seq and stores those in dict.

    filename: str, name of argonaut formatted input file
    return: nested dictionary with {accession_number:{organism_name:DNA-seq}}
    """

    lines=(open(filename))
    record_bool = False
    records=[]
    record = ""
    for line in lines:
        line=line.strip()
        line = line.replace('\n', '')
        line = line.replace('\t', '')
        if line.startswith("BEGIN IONS"):
            #line = line.replace("ACCESSION   ", "")
            record_bool=True
        elif line.startswith("END IONS"):
            record_bool=False
            records.append(record)
            record=""
        # elif line.startswith("ORGANISM"):
        #     organism_dict = {}
        #     line = line.replace("ORGANISM  ", "")
        #     key = list(gb_dict)[-1]
        #     organism_dict[line] = ""
        #     gb_dict[key] = organism_dict
        # elif "ORIGIN" in line:
        #     origin=True
        # elif "//" in line:
        #     origin=False
        if record_bool:
            record+=line
    #print(records)

    # dict={}
    for record in records:
        key = re.search(r'190.119(.*)520', record)
        if key is not None:
            print(record)
        # key=re.search(r'NAME=(.*)',record)
        # smiles=re.search(r'SMILES=(.*)',record)
        # mass=re.search(r'PEPMASS=(.*)',record)
        # if key is not None:
        #     dict[key.group(1)]=[]
        # if smiles is not None:
        #     last_key = list(dict)[-1]
        #     dict[last_key].append(smiles.group(1))
        # if mass is not None:
        #     last_key = list(dict)[-1]
        #     dict[last_key].append(mass.group(1))
        #INCHI=
    # print(dict)
            # line=line.replace(" ", "")
            # line=''.join(filter(lambda ch: not ch.isdigit(), line))
            # #https://www.studytonight.com/python-howtos/remove-numbers-from-string-in-python
            # line = line.replace("ORIGIN", "")
            # key = list(gb_dict)[-1]
            # dict=gb_dict[key]
            # last_key=list(dict)[-1]
            # gb_dict[key][last_key]+=line
    return None

def read_json(json_file):
    """

    :param json_file:
    :return:
    """
    with open(json_file, 'r') as f:
        dict=json.load(f)
    df=json_normalize(dict)
    df.set_index("spectrum_id",inplace=True, drop=True)
    #print(df.loc["CCMSLIB00004678842"])
    #print(df.loc["CCMSLIB00000072521", "peaks_json"])
    return df

def new_dataframe(df_massql_matches,df_json):
    """
    Makes a dataframe with scan precmz smiles and
    :return:
    """
    # this is also an interesting column for df_json but I think its always the same as the smiles "InChIKey_smiles"
    df=pd.merge(df_massql_matches["precmz"],df_json[["Precursor_MZ","Smiles", "peaks_json"]],left_index=True, right_index=True)
    # for some matches there are no smiles so remove those from the dataframe
    df.drop(df.index[df['Smiles'] == 'N/A'], inplace=True)
    df.drop(df.index[df['Smiles'] == ' '], inplace=True)
    print(df)
    return df

def make_spectrum_file_for_id(df_json, spectrum_id, path_to_store_spectrum_files):
    """Takes a nested sorted list and outputs a tab delimited file

    alignment_list: nested list with families and alignment lenghts
    return: tab delimited text file with the contents of each sub list on a line
    """
    list_of_lists=ast.literal_eval(df_json.loc[spectrum_id, "peaks_json"])
    file_path = Path(r"{0}/spectrum_file_{1}.txt".format(path_to_store_spectrum_files,spectrum_id))
    # if os.path.exists(file_path):
    #     print(os.path.abspath(file_path))
    #     return os.path.abspath(file_path)
    spectrum_file=open(file_path, "w")
    for i in range(3):
        spectrum_file.write("energy{0}\n".format(i))
        for sub_list in list_of_lists:
            if sub_list[1]>600: #dit is ff een tussen oplossing om een file te krijgen waar je iets mee kan!
                spectrum_file.write("{0} {1}".format(sub_list[0], sub_list[1]))
                spectrum_file.write("\n")
    spectrum_file.close()
    print(os.path.abspath(file_path))

def make_spectrum_file_for_id_matchms(gnps_pickled_lib, spectrum_id, path_to_store_spectrum_files):
    #object=list(matchms.importing.load_scores.scores_from_pickle(gnps_pickled_lib))
    obj = pd.read_pickle(gnps_pickled_lib)
    # with open(gnps_pickled_lib, 'r') as f:
    #     spectra=pickletools.(f)
    # spectra=gnps_pickled_libS
    # print(spectra.peaks.mz[0])
    for spectrum in obj:
        if spectrum.get("spectrumid") == spectrum_id:
            file_path = Path(r"{0}/spectrum_file_{1}_new.txt".format(path_to_store_spectrum_files, spectrum_id))
            print(os.path.abspath(file_path))
            spectrum_file=open(file_path, "w")
            for i in range(3):
                spectrum_file.write("energy{0}\n".format(i))
                for fragment in range(len(spectrum.peaks.mz)):
                    spectrum_file.write("{0} {1}".format(spectrum.peaks.mz[fragment], spectrum.peaks.intensities[fragment]))
                    spectrum_file.write("\n")
            spectrum_file.close()
    #object=list(matchms.importing.load_from_pickle)
    #print(object)
    return None

def annotate_peaks(spectrum_file, smiles):
    """
    Annotates the MS2 peaks given a smiles and a fragmentation spectrum
    :return:
    """
    #print(df.loc["CCMSLIB00006126912"])
    return None

def annotate_peaks(spectrum_file_name, smiles, identifier, abs_mass_tol=0.01):
    """
    Annotates the MS2 peaks given a smiles and a fragmentation spectrum
    :return:
    """
    # print(df.loc["CCMSLIB00004678842"])
    # cfm-annotate.exe <smiles_or_inchi> <spectrum_file> **<id>** <ppm_mass_tol> <abs_mass_tol> 0.01<param_file> <config_file> <output_file
    file_path_out = Path(r"/lustre/BIF/nobackup/seele006/cfm_annotation_out_{0}".format(identifier))
    out_fn = "cfm_annotation_out_{0}".format(identifier)
    if os.path.exists(out_fn):
        return out_fn
    cmd = 'cfm-annotate \'{0}\' {1} {2} -abs_mass_tol {3} -output_file {4}' \
            .format(smiles, spectrum_file_name, identifier, abs_mass_tol, out_fn)
    e = subprocess.check_call(cmd, shell=True)
    print("EXIT STATUS AND TYPE", e, type(e))
    print("hi")
    return out_fn

def visualize_mol(smiles: str) -> None:
    """
    Takes a smiles as input and outputs an PIL<PNG> image

    :param smiles: str, a smiles string that corresponds to a molecular structure
    :return: PIL<PNG> image, image of structure corresponding with the smiles

    function adapted from: https://pchanda.github.io/See-substructure-in-molecule/
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    Chem.Kekulize(mol)
    img = MolToImage(mol, size=(200, 200), fitImage=True)
    return img

def look_for(cfm_file,motif,fragments):
    """
    # First parse the cfm_file because you only want the middle part
    #Then for each fragment of the motif you want you want to select the smiles
    Selection of the annotated peaks that are relevant for the mass2motif
    #TODO: think about how you will do this for neutral loss...
    :return:
    """
    return None

def select_the_most_likely_option_for_one_motif():
    #look through all the possible smiles for a fragment and select the one that is there most often
    # return the motif with the with fragment, probability and most likely smiles
    return None

def delete_files(path_to_store_spectrum_files):
    for file in os.scandir(path_to_store_spectrum_files):
        os.remove(file.path)
    return None

def main():
    """Main function of this module"""
    path_to_pickle_file = argv[2]
    filename=argv[1]
    path_to_store_spectrum_files = argv[3]
    path_to_store_json_file=argv[4]
    make_json_file(path_to_pickle_file, path_to_store_json_file)
    # step 0: parse input line
    lines = (open(filename))
    for line in lines:
        line = line.strip()
        line = line.replace('\n', '')
        motif, fragments, query=parse_line_with_motifs_and_querries(line)
        query = ("QUERY scaninfo(MS2DATA) WHERE POLARITY = Positive AND MS2PROD = 68.0275:TOLERANCEMZ=0.01 AND MS2PROD = 85.0250:TOLERANCEMZ=0.01 AND MS2PROD = 97.0250:TOLERANCEMZ=0.01")
        # step 1: parse json file
        df_json=read_json(path_to_store_json_file)
        # step 2: search query in json file with MassQL
        df_massql_matches=try_massql(query, path_to_store_json_file)
        # step 3: get the smiles for every match of MassQL
        df_matches_and_smiles=new_dataframe(df_massql_matches,df_json)
        # step 4: print a spectrum file for a match
        #for identifier in list(index_smiles)
        #list_of_lists = ast.literal_eval(df_json.loc[identifier, "peaks_json"])
        #make_spectrum_file_for_id(list_of_lists, identifier)
        identifier="CCMSLIB00000001645"
        spectrum_file_name=make_spectrum_file_for_id(df_json, identifier, path_to_store_spectrum_files)
        #print(df_matches_and_smiles.loc[identifier, "Smiles"])
        #make_spectrum_file_for_id_matchms(path_to_pickle_file, identifier, path_to_store_spectrum_files)
        # make a huge list for each of the motifs containing the possible smiles per fragments
        #annotate_peaks(spectrum_file_name, smiles)
        #Make PDF
        # pdf = FPDF()
        # pdf.add_page()
        # pdf.set_font("helvetica", size=10)
        # pdf.image(visualize_mol("N[C@@H](CCCCNC(N)=O)C(O)=O"))
        # pdf.output("output.pdf")


if __name__ == "__main__":
    main()