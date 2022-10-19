#!/usr/bin/env massql
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
import json
from typing import List
import numpy
from matchms import Spectrum
import csv
from matchms.importing import load_from_mgf
from matchms.exporting import save_as_mgf

# functions

def parse_line_with_motifs_and_querries(line):
    motif=re.search(r'(.*)    (.*)    (.*)', line).group(1)
    fragments=re.search(r'(.*)    (.*)    (.*)', line).group(2)
    query=re.search(r'(.*)    (.*)    (.*)', line).group(3)
    return motif,fragments,query

def try_massql(query, json_file):
    print(query)
    # this msql_engine only takes json files
    df=msql_engine.process_query(query,json_file)
    print(df)
    df.rename(columns={'scan': 'spectrum_id'}, inplace=True)
    df.set_index("spectrum_id", inplace=True, drop=True)
    return df

def read_json(json_file):
    """

    :param json_file:
    :return:
    """
    with open(json_file, 'r') as f:
        dict=json.load(f)
    df=json_normalize(dict)
    #print(df.columns)
    df.set_index("spectrum_id",inplace=True, drop=True)
    #print(df.loc["CCMSLIB00004678842"])
    #print(df.loc["CCMSLIB00000072521", "peaks_json"])
    print(df)
    return df

def new_dataframe(df_massql_matches,df_json):
    """
    Makes a dataframe with scan precmz smiles and
    :return:
    """
    # this is also an interesting column for df_json but I think its always the same as the smiles "InChIKey_smiles"
    df=pd.merge(df_massql_matches["precmz"],df_json[["Precursor_MZ","Smiles", "peaks_json"]],left_index=True, right_index=True)
    # for pickle file smiles instead of Smiles
    # for some matches there are no smiles so remove those from the dataframe
    df.drop(df.index[df['Smiles'] == 'N/A'], inplace=True)
    df.drop(df.index[df['Smiles'] == ' '], inplace=True)
    # for index, row in df.iterrows():
    #     print(df.at[index,"Smiles"])
    print(df)
    return df

def make_spectrum_file_for_id2(motif, df_json, spectrum_id, path_to_store_spectrum_files, mgf_file):
    """Takes a nested sorted list and outputs a tab delimited file

    alignment_list: nested list with families and alignment lenghts
    return: tab delimited text file with the contents of each sub list on a line
    """
    #list_of_lists=ast.literal_eval(df_json.loc[spectrum_id, "peaks_json"])
    HMDB_id = re.search(r'(HMDB:)(HMDB\d*)(-\d*)(.*)', df_json.loc[spectrum_id, "Compound_Name"]).group(2)
    # in the current (10-2022) structures database of HMDB a longer identifier is used, so the older identifier from GNPS
    # needs to be adjusted
    adj_HMDB_id = HMDB_id[:3] + '00' + HMDB_id[3:]
    file_path = Path(fr"{path_to_store_spectrum_files}/spectrum_{motif}_{adj_HMDB_id}.txt")

    with open(mgf_file, 'r') as spectra_file:
        spectra_from_file = list(load_from_mgf(spectra_file))
        for spectrum in spectra_from_file:
            print(spectrum.get("spectrumid"))
            if spectrum.get("spectrumid") == spectrum_id:
                print(spectrum.get("spectrum_id"))
                save_as_mgf(spectrum, file_path)

    # if os.path.exists(file_path):
    #     print(os.path.abspath(file_path))
    #     return os.path.abspath(file_path)
    # spectrum_file=open(file_path, "w")
    # for i in range(3):
    #     spectrum_file.write("energy{0}\n".format(i))
    #     for sub_list in list_of_lists:
    #         #if sub_list[1]>600: #dit is ff een tussen oplossing om een file te krijgen waar je iets mee kan!
    #         spectrum_file.write("{0} {1}".format(sub_list[0], sub_list[1]))
    #         spectrum_file.write("\n")
    # spectrum_file.close()
    print(os.path.abspath(file_path))

def mgf_to_spectra(mgf_file):
    with open(mgf_file, 'r') as spectra_file:
        spectra_from_file = list(load_from_mgf(spectra_file))
        for spectrum in spectra_from_file:
            if spectrum.get("spectrum_id") == spectrum_id:

                save_as_mgf(spectrum, filename)


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
            #if sub_list[1]>600: #dit is ff een tussen oplossing om een file te krijgen waar je iets mee kan!
            spectrum_file.write("{0} {1}".format(sub_list[0], sub_list[1]))
            spectrum_file.write("\n")
    spectrum_file.close()
    print(os.path.abspath(file_path))

def delete_files(path_to_store_spectrum_files):
    for file in os.scandir(path_to_store_spectrum_files):
        os.remove(file.path)
    return None

def main():
    """Main function of this module"""
    #path_to_pickle_file = argv[2]
    path_to_json_file= argv[2]
    filename=argv[1]
    path_to_store_spectrum_files = argv[3]
    mgf_file= argv[4]
    #path_to_json_file=argv[4]
    #make_json_file(path_to_pickle_file, path_to_json_file)
    # step 0: parse input line
    lines = (open(filename))
    #save_json_as_csv(path_to_json_file)
    for line in lines:
        line = line.strip()
        line = line.replace('\n', '')
        motif, fragments, query=parse_line_with_motifs_and_querries(line)
        #query = ("QUERY scaninfo(MS2DATA) WHERE POLARITY = Positive AND MS2PROD = 85.0250:TOLERANCEMZ=0.01:INTENSITYMATCH=Y:INTENSITYMATCHREFERENCE AND MS2PROD = 68.0275:TOLERANCEMZ=0.01:INTENSITYMATCH=Y*0.186:INTENSITYMATCHPERCENT=99 AND MS2PROD = 97.0250:TOLERANCEMZ=0.01:INTENSITYMATCH=Y*0.156:INTENSITYMATCHPERCENT=99")
        #query = ("QUERY scaninfo(MS2DATA) WHERE POLARITY = Positive AND MS2NL = 176.0350:TOLERANCEMZ=0.01:INTENSITYMATCH=Y:INTENSITYMATCHREFERENCE AND MS2PROD = 126.0550:TOLERANCEMZ=0.01:INTENSITYMATCH=Y*0.089:INTENSITYMATCHPERCENT=99 AND MS2PROD = 127.0375:TOLERANCEMZ=0.01:INTENSITYMATCH=Y*0.082:INTENSITYMATCHPERCENT=99")
        #query = ("QUERY scaninfo(MS2DATA) WHERE POLARITY = Positive AND MS2NL = 46.0050:TOLERANCEMZ=0.005") #motif gnps_motif_38.m2m
        # step 1: parse json file
        df_json=read_json(path_to_json_file)
        # step 2: search query in json file with MassQL
        df_massql_matches=try_massql(query, path_to_json_file)
        # step 3: get the smiles for every match of MassQL
        df_matches_and_smiles=new_dataframe(df_massql_matches,df_json)

        # step 4: print a spectrum file for a match
        #for identifier in list(index_smiles)
        #list_of_lists = ast.literal_eval(df_json.loc[identifier, "peaks_json"])
        #make_spectrum_file_for_id(list_of_lists, identifier)

        #for identifier, row in df_matches_and_smiles.iterrows():
        #make a spectrum file
        identifier="CCMSLIB00000426038" #result from HMDB with Motif_38
        spectrum_file_name=make_spectrum_file_for_id2(motif, df_json, identifier, path_to_store_spectrum_files, mgf_file)
        print(df_matches_and_smiles.loc[identifier, "Smiles"])
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