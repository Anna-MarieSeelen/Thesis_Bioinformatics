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
from pandas import json_normalize
import pandas as pd
from pathlib import Path
import re
import os
import json
from matchms.importing import load_from_mgf
from matchms.exporting import save_as_mgf

# functions

def parse_input(filename):
    """Parses mgf into strings, where each string is a spectrum and stores those in dict with the spectrum_id as key.

    filename: str, name of mgf formatted input file
    return: dictionary with {spectrum_id:record} where each record is a string containing the mgf-style spectrum of one
    compound
    """

    lines=(open(filename))
    record_bool = False
    records=[]
    record = ""
    for line in lines:
        if line.startswith("BEGIN IONS"):
            record_bool=True
        if record_bool:
            record+=line
        if line.startswith("END IONS"):
            records.append(record)
            record_bool = False
            record=""

    dict_with_mgf_spectra={}
    for record in records:
        record=record.strip()
        # use look for the spectrumid in the string and use it as a key for the dict
        key = re.search(r'SPECTRUMID=(.*)', record).group(1)
        if key is not None:
            dict_with_mgf_spectra[key] = record
    return dict_with_mgf_spectra

def parse_line_with_motifs_and_querries(line: str) -> tuple:
    """
    Takes a line in a tab separated file and separates it into motif, feature list and massql query

    :param line: str, line in a tab separated file containing the motif, list of features and the massql query
    :return: returns the motif, feature_list and the massql query in a tuple
    """
    motif=re.search(r'(.*)    (.*)    (.*)', line).group(1)
    fragments=re.search(r'(.*)    (.*)    (.*)', line).group(2)
    query=re.search(r'(.*)    (.*)    (.*)', line).group(3)
    return motif,fragments,query

def try_massql(query, json_file):
    """
    Uses MassQL to retrieve the identifiers of the spectra in the json file that contain the characteristics of the query

    :param query: str, MassQL query made based on a mass2motif
    :param json_file: str, the path to a file with MS/MS spectra retrieved from GNPS in json format
    :return: returns a Pandas dataframe with the spectrum ids of the spectra as the index which contain the
    characteristics of the query.
    """
    print(query)
    # this msql_engine only takes json files
    df=msql_engine.process_query(query,json_file)
    print(df)
    df.rename(columns={'scan': 'spectrum_id'}, inplace=True)
    df.set_index("spectrum_id", inplace=True, drop=True)
    return df

def read_json(json_file):
    """
    Puts a json file with MS/MS spectra retrieved in a Pandas dataframe.

    :param json_file: str, the path to a file with MS/MS spectra retrieved from GNPS in json format
    :return: pandas Dataframe with one row for each spectrum and the information belonging to the spectrum in separate
    columns.
    """
    with open(json_file, 'r') as f:
        dict=json.load(f)
    df=json_normalize(dict)
    #print(df.columns)
    df.set_index("spectrum_id",inplace=True, drop=True)
    print(df)
    return df

def make_spectrum_file_for_id2(motif, df_json, spectrum_id, path_to_store_spectrum_files, dict_with_mgf_spectra):
    """Takes a nested sorted list and outputs a tab delimited file

    alignment_list: nested list with families and alignment lenghts
    return: tab delimited text file with the contents of each sub list on a line
    """
    # look for the corresponding HMDB identifier using the spectrum_id in the dataframe of the json file,
    # because MAGMa works with that identifier.
    #TODO: delete the json read function and don't use the dataframe anymore, just read this from the dict string, it will be easy.
    #spectrum=dict_with_mgf_spectra[spectrum_id]
    HMDB_id = re.search(r'(HMDB:)(HMDB\d*)(-\d*)(.*)', df_json.loc[spectrum_id, "Compound_Name"]).group(2)
    # HMDB_id = re.search(r'(.*)(HMDB:)(HMDB\d*)(-\d*)(.*)', spectrum).group(3)
    # in the current (10-2022) structures database of HMDB a longer identifier is used, so the older identifier from GNPS
    # needs to be adjusted
    adj_HMDB_id = HMDB_id[:4] + '00' + HMDB_id[4:]
    file_path = Path(fr"{path_to_store_spectrum_files}/spectrum_{motif}_{adj_HMDB_id}.txt")
    if os.path.exists(file_path):
         return os.path.abspath(file_path)
    else:
        spectrum_file = open(file_path, "w")
        spectrum_file.write(dict_with_mgf_spectra[spectrum_id])
        spectrum_file.close()
        return os.path.abspath(file_path)

def write_path_to_file(path_to_files_with_motifs, path_to_spectrum_file):
    # the path to the spectrum files that are associated with all motifs according to MassQL will be stored in this document
    path, filename = os.path.split(path_to_files_with_motifs)
    file_path_out = Path(fr"{path}/paths_to_spectrum_files_from_MassQL.txt")
    if os.path.exists(file_path_out):
        output_file = open(file_path_out, "a")
        output_file.write(f"{path_to_spectrum_file}")
        output_file.write("\n")
        output_file.close()
    else:
        output_file = open(file_path_out, "w")
        output_file.write(f"{path_to_spectrum_file}")
        output_file.write("\n")
        output_file.close()
    return None

def delete_files(path_to_store_spectrum_files):
    for file in os.scandir(path_to_store_spectrum_files):
        os.remove(file.path)
    return None

def main():
    """Main function of this module"""
    #path_to_pickle_file = argv[2]
    path_to_json_file= argv[2]
    path_to_files_with_motifs=argv[1]
    path_to_store_spectrum_files = argv[3]
    mgf_file= argv[4]
    dict_with_mgf_spectra=parse_input(mgf_file)
    #path_to_json_file=argv[4]
    #make_json_file(path_to_pickle_file, path_to_json_file)
    # step 0: parse input line
    lines = (open(path_to_files_with_motifs))
    #save_json_as_csv(path_to_json_file)
    for line in lines:
        line = line.strip()
        line = line.replace('\n', '')
        motif, fragments, query=parse_line_with_motifs_and_querries(line)
        #query = ("QUERY scaninfo(MS2DATA) WHERE POLARITY = Positive AND MS2PROD = 85.0250:TOLERANCEMZ=0.01:INTENSITYMATCH=Y:INTENSITYMATCHREFERENCE AND MS2PROD = 68.0275:TOLERANCEMZ=0.01:INTENSITYMATCH=Y*0.186:INTENSITYMATCHPERCENT=99 AND MS2PROD = 97.0250:TOLERANCEMZ=0.01:INTENSITYMATCH=Y*0.156:INTENSITYMATCHPERCENT=99")
        #query = ("QUERY scaninfo(MS2DATA) WHERE POLARITY = Positive AND MS2NL = 176.0350:TOLERANCEMZ=0.01:INTENSITYMATCH=Y:INTENSITYMATCHREFERENCE AND MS2PROD = 126.0550:TOLERANCEMZ=0.01:INTENSITYMATCH=Y*0.089:INTENSITYMATCHPERCENT=99 AND MS2PROD = 127.0375:TOLERANCEMZ=0.01:INTENSITYMATCH=Y*0.082:INTENSITYMATCHPERCENT=99")
        #query = ("QUERY scaninfo(MS2DATA) WHERE POLARITY = Positive AND MS2NL = 46.0050:TOLERANCEMZ=0.005") #motif gnps_motif_38.m2m
        # step 2: search query in json file with MassQL
        df_massql_matches=try_massql(query, path_to_json_file)
        # step 1: to make it easy to get the smiles belonging to the matches and the HMDB identifiers by putting the
        # json file in a Pandas Dataframe
        df_json=read_json(path_to_json_file)
        # step 3: get the smiles for every match of MassQL
        df_matches_and_smiles=new_dataframe(df_massql_matches,df_json)

        for identifier, row in df_matches_and_smiles.iterrows():
            #make a spectrum file
            #identifier="CCMSLIB00000426038" #result from HMDB with Motif_38
            path_to_spectrum_file=make_spectrum_file_for_id2(motif, df_json, identifier, path_to_store_spectrum_files, dict_with_mgf_spectra)
            # add to list of spectrum file
            print(df_matches_and_smiles.loc[identifier, "Smiles"])
            write_path_to_file(path_to_files_with_motifs, path_to_spectrum_file)


if __name__ == "__main__":
    main()