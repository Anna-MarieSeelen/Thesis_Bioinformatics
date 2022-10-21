#!/usr/bin/env massql
"""
Author: Anna-Marie Seelen
Studentnumber:1008970
Description: Takes a file with MassQL queries, searches for spectra in json file that contain the queries and writes
these spectra to individual files in mgf-style.
Usage: python3 *massql_test.py *path_to_file_with_motifs *path_to_HMDB_json_file* *path_to_store_spectrum_files*
*path_to_HMDB_mgf_file*
    path_to_file_with_motifs: a tab separated file with a selected motif, feature list and massql query on each line
    (output from make_pdf_with_smiles.py)
    path_to_HMDB_json_file: the path to the file with HMDB MS/MS spectra retrieved from GNPS in json format
    path_to_store_spectrum_files: folder where all the mgf-formatted text files with spectra will be stored that contain
    a selected motif (determined by MassQL)
    path_to_HMDB_mgf_file: the path to the file with mgf-formatted HMDB MS/MS spectra from GNPS
"""

# import statements
from massql import msql_engine
from sys import argv
from pathlib import Path
import re
import os
import pandas as pd

# functions
def parse_input(HMDB_mgf_file: str) -> dict:
    """Parses mgf into strings, where each string is a spectrum and stores those in dict with the spectrum_id as key.

    :param HMDB_mgf_file: str, name of mgf formatted file containing MS/MS spectra from GNPS
    :return: dictionary with {spectrum_id:record} where each record is a string containing the mgf-style accession of one
    compound
    """

    lines_HMDB_mgf_file=open(HMDB_mgf_file)
    spectrum_record_bool = False
    mgf_spectrum_records=[]
    mgf_spectrum_record = ""
    for line in lines_HMDB_mgf_file:
        if line.startswith("BEGIN IONS"):
            spectrum_record_bool=True
        if spectrum_record_bool:
            mgf_spectrum_record+=line
        if line.startswith("END IONS"):
            mgf_spectrum_records.append(mgf_spectrum_record)
            spectrum_record_bool = False
            mgf_spectrum_record=""

    dict_with_mgf_spectra={}
    for mgf_spectrum_record in mgf_spectrum_records:
        mgf_spectrum_record=mgf_spectrum_record.strip()
        # look for the spectrumid in the string and use it as a key for the dict
        key = re.search(r'SPECTRUMID=(.*)', mgf_spectrum_record).group(1)
        if key is not None:
            dict_with_mgf_spectra[key] = mgf_spectrum_record
    return dict_with_mgf_spectra

def parse_line_with_motif_and_query(line: str) -> tuple:
    """
    Takes a line in a tab separated file and separates it into motif, feature list and massql query

    :param line: str, line in a tab separated file containing the motif, list of features and the massql query
    :return: returns the motif, feature_list and the massql query in a tuple
    """
    motif=re.search(r'(.*)    (.*)    (.*)', line).group(1)
    features=re.search(r'(.*)    (.*)    (.*)', line).group(2)
    massql_query=re.search(r'(.*)    (.*)    (.*)', line).group(3)
    return motif,features,massql_query

class Suppress:
    """
    Class to supress the print statement in the msql_engine code
    """
    def __init__(self, *, suppress_stdout=False, suppress_stderr=False):
        self.suppress_stdout = suppress_stdout
        self.suppress_stderr = suppress_stderr
        self.original_stdout = None
        self.original_stderr = None

    def __enter__(self):
        import sys, os
        devnull = open(os.devnull, "w")

        # Suppress streams
        if self.suppress_stdout:
            self.original_stdout = sys.stdout
            sys.stdout = devnull

        if self.suppress_stderr:
            self.original_stderr = sys.stderr
            sys.stderr = devnull

    def __exit__(self, *args, **kwargs):
        import sys
        # Restore streams
        if self.suppress_stdout:
            sys.stdout = self.original_stdout

        if self.suppress_stderr:
            sys.stderr = self.original_stderr

def search_motif_with_massql(massql_query: str, path_to_HMDB_json_file: str) -> pd.DataFrame:
    """
    Uses MassQL to retrieve the identifiers of the spectra in the json file that contain the query

    :param massql_query: str, MassQL query made based on a Mass2Motif
    :param path_to_HMDB_json_file: str, the path to a file with MS/MS spectra retrieved from GNPS in json format
    :return: returns a Pandas dataframe with the spectrum ids of the spectra as the index which contain the
    characteristics of the query.
    """
    # msql_engine only takes json files
    # the class Suppress suppresses the output to the command line that the msql_engine is giving through a print
    # function.
    with Suppress(suppress_stderr=True, suppress_stdout=True):
        df_massql_res=msql_engine.process_query(massql_query,path_to_HMDB_json_file)
    df_massql_res.rename(columns={'scan': 'spectrum_id'}, inplace=True)
    df_massql_res.set_index("spectrum_id", inplace=True, drop=True)
    return df_massql_res

def make_mgf_txt_file_for_spectrum(motif: str, spectrum_id: str, path_to_store_spectrum_files: str,
                                       dict_with_mgf_spectra: dict) -> str:
    """
    Writes a spectrum containing the motif to a text file in mgf format using the identifier given by MassQL

    :param motif: str, the Mass2Motif for which the MassQL query was made and for which the spectrum was found, which
    contains the motif
    :param spectrum_id: str, the spectrum_id of GNPS spectra files formatted like this: CCMSLIB00000425029
    :param path_to_store_spectrum_files: str, folder where all the mgf-formatted text files with spectra will be stored.
    :param dict_with_mgf_spectra: dict, with {spectrum_id:record} where each record is a string containing the mgf-style
    spectrum of one compound
    :return: path of the spectrum file with the HMDB identifier of the spectrum and the motif it is supposed to contain
    in the file name.
    """
    # look for the corresponding HMDB identifier using the spectrum_id in dict of the mgf file,
    # because MAGMa works with the HMDB identifier.
    spectrum_record_mgf=dict_with_mgf_spectra[spectrum_id]
    HMDB_id_of_spectrum = re.search(r'(.*)(HMDB:)(HMDB\d*)(-\d*)(.*)', spectrum_record_mgf).group(3)
    # in the current (10-2022) structures database of HMDB which is used for MAGMa a longer identifier is used, so the
    # older identifier from GNPS needs to be adjusted
    adj_HMDB_id = HMDB_id_of_spectrum[:4] + '00' + HMDB_id_of_spectrum[4:]
    # create the name for the spectrum file with the HMDB identifier and the motif that the spectrum contains according
    # to MassQL
    path_to_spectrum_file = Path(fr"{path_to_store_spectrum_files}/spectrum_{motif}_{adj_HMDB_id}.txt")
    # write the mgf spectrum saved in the dict only to a text file if the text file doesn't exist already
    if os.path.exists(path_to_spectrum_file):
         return os.path.abspath(path_to_spectrum_file)
    else:
        spectrum_file = open(path_to_spectrum_file, "w")
        spectrum_file.write(dict_with_mgf_spectra[spectrum_id])
        spectrum_file.close()
        return os.path.abspath(path_to_spectrum_file)

def write_path_to_file(path_to_file_with_motifs: str, path_to_spectrum_file: str) -> None:
    """
    Makes or appends a path of a spectrum file to a text file with paths of the spectrum files which contain Mass2Motifs

    :param path_to_file_with_motifs: str, the path where the file with the selected motifs is stored
    :param path_to_spectrum_file: str, path of the spectrum file with the HMDB identifier of the spectrum and the motif
    it is supposed to contain in the file name.
    :return:
    """
    # the paths to the spectrum files that are associated with all motifs according to MassQL will be stored in this
    # document. This document will be used by the MAGMa.py script.
    path, filename = os.path.split(path_to_file_with_motifs)

    path_to_txt_file_with_paths = Path(fr"{path}/paths_to_spectrum_files_from_MassQL.txt")
    # if the path to the spectrum names file already exists you should append to this file
    if os.path.exists(path_to_txt_file_with_paths):
        # however, you should only append a path that is not yet in the file!
        with open(path_to_txt_file_with_paths, 'r') as file:
            if path_to_spectrum_file in file.read():
                file.close()
                return None
            else:
                output_file = open(path_to_txt_file_with_paths, "a")
                output_file.write(f"{path_to_spectrum_file}")
                output_file.write("\n")
                output_file.close()
    # if the path to the spectrum names file doesn't exist make one and write the path to the spectrum file in it.
    else:
        output_file = open(path_to_txt_file_with_paths, "w")
        output_file.write(f"{path_to_spectrum_file}")
        output_file.write("\n")
        output_file.close()
    return None

def main():
    """Main function of this module"""
    path_to_file_with_motifs = argv[1]
    path_to_HMDB_json_file = argv[2]
    path_to_store_spectrum_files = argv[3]
    path_to_HMDB_mgf_file = argv[4]
    # step 1: parse the mgf file with spectra and put into dictionary with spectrum id's as keys.
    dict_with_mgf_spectra=parse_input(path_to_HMDB_mgf_file)
    # step 2: parse the lines in the file where all the selected motifs and their corresponding massql queries are
    # listed.
    lines_motif_file = open(path_to_file_with_motifs)
    for line in lines_motif_file:
        line = line.strip()
        line = line.replace('\n', '')
        # step 3: retrieve the motif and the massql query one by one
        motif, features, massql_query=parse_line_with_motif_and_query(line)
        # step 4: search for spectra with the motif in json file with MassQL
        df_massql_matches=search_motif_with_massql(massql_query, path_to_HMDB_json_file)
        for identifier, row in df_massql_matches.iterrows():
            #identifier="CCMSLIB00000426038" #result from HMDB with Motif_38
            # step 5: make a spectrum file in mgf format for every identified match
            path_to_spectrum_file = make_mgf_txt_file_for_spectrum(motif, identifier, path_to_store_spectrum_files,
                                                               dict_with_mgf_spectra)
            # step 6: add the new spectrum file to a document with a list of spectrum files to be used by MAGMa for
            # annotation.
            write_path_to_file(path_to_file_with_motifs, path_to_spectrum_file)

if __name__ == "__main__":
    main()