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
from pathlib import Path
import re
import os

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
    # this msql_engine only takes json files
    df=msql_engine.process_query(query,json_file)
    df.rename(columns={'scan': 'spectrum_id'}, inplace=True)
    df.set_index("spectrum_id", inplace=True, drop=True)
    return df

def make_spectrum_file_for_id2(motif, spectrum_id, path_to_store_spectrum_files, dict_with_mgf_spectra):
    """
    Writes a spectrum containing the motif to a text file in mgf format using the identifier given by MassQL

    :param motif: str, the Mass2Motif for which the MassQL query was made and for which the spectrum was found, which
    contains the motif
    :param spectrum_id: str, The spectrum_id of GNPS CCMSLIB00000425029
    :param path_to_store_spectrum_files: str, folder where all the mgf-formatted text files with spectra will be stored.
    :param dict_with_mgf_spectra: dict, with {spectrum_id:record} where each record is a string containing the mgf-style
    spectrum of one compound
    :return: path of the spectrum file with the HMDB identifier of the spectrum and the motif it is supposed to contain
    in the file name.
    """
    # look for the corresponding HMDB identifier using the spectrum_id in dict of the mgf file,
    # because MAGMa works with the HMDB identifier.
    spectrum=dict_with_mgf_spectra[spectrum_id]
    HMDB_id = re.search(r'(.*)(HMDB:)(HMDB\d*)(-\d*)(.*)', spectrum).group(3)
    # in the current (10-2022) structures database of HMDB which is used for MAGMa a longer identifier is used, so the
    # older identifier from GNPS needs to be adjusted
    adj_HMDB_id = HMDB_id[:4] + '00' + HMDB_id[4:]
    # create the name for the spectrum file
    file_path = Path(fr"{path_to_store_spectrum_files}/spectrum_{motif}_{adj_HMDB_id}.txt")
    # write the mgf spectrum saved in the dict only to a text file if the text file doesn't exist already
    if os.path.exists(file_path):
         return os.path.abspath(file_path)
    else:
        spectrum_file = open(file_path, "w")
        spectrum_file.write(dict_with_mgf_spectra[spectrum_id])
        spectrum_file.close()
        return os.path.abspath(file_path)

def write_path_to_file(path_to_files_with_motifs, path_to_spectrum_file):
    """
    Makes or appends a path of a spectrum file to a text file with paths of the spectrum files which contain Mass2motifs

    :param path_to_files_with_motifs: str, the path where the file with the selected motifs is stored
    :param path_to_spectrum_file: str, path of the spectrum file with the HMDB identifier of the spectrum and the motif
    it is supposed to contain in the file name.
    :return:
    """
    # the path to the spectrum files that are associated with all motifs according to MassQL will be stored in this
    # document. This document will be used by the MAGMa.py script.
    path, filename = os.path.split(path_to_files_with_motifs)
    file_path_out = Path(fr"{path}/paths_to_spectrum_files_from_MassQL.txt")
    # if the path to the spectrum names file already exists you should append to this file
    if os.path.exists(file_path_out):
        # however, you should only append a path that is not yet in the file!
        with open(file_path_out, 'r') as file:
            if path_to_spectrum_file in file.read():
                file.close()
                return None
            else:
                output_file = open(file_path_out, "a")
                output_file.write(f"{path_to_spectrum_file}")
                output_file.write("\n")
                output_file.close()
    # if the path to the spectrum names file doesn't exist make one and write the path to the spectrum file in it.
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
    path_to_json_file= argv[2]
    path_to_files_with_motifs=argv[1]
    path_to_store_spectrum_files = argv[3]
    mgf_file= argv[4]
    dict_with_mgf_spectra=parse_input(mgf_file)
    # step 0: parse input line
    lines = (open(path_to_files_with_motifs))
    for line in lines:
        line = line.strip()
        line = line.replace('\n', '')
        motif, fragments, query=parse_line_with_motifs_and_querries(line)
        # step 2: search query in json file with MassQL
        df_massql_matches=try_massql(query, path_to_json_file)
        print(df_massql_matches)
        # step 3: get the smiles for every match of MassQL
        #df_matches_and_smiles=new_dataframe(df_massql_matches,df_json)
        # for identifier (so actually every spectrum) in the HMDB file that contains the motif/MassQL query:
        for identifier, row in df_massql_matches.iterrows():
            #identifier="CCMSLIB00000426038" #result from HMDB with Motif_38
            # make a spectrum file in mgf format
            path_to_spectrum_file=make_spectrum_file_for_id2(motif, identifier, path_to_store_spectrum_files, dict_with_mgf_spectra)
            # add the new spectrum file to a document with a list of spectrum files
            write_path_to_file(path_to_files_with_motifs, path_to_spectrum_file)

if __name__ == "__main__":
    main()