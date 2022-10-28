#!/usr/bin/env python3
"""
Author: Anna-Marie Seelen
Studentnumber:1008970
Description: Takes a file with MassQL queries, searches for spectra in json file that contain the queries and writes
these spectra to individual files in mgf-style.
Usage: python3 *massql_search_spectra_with_motif.py *path_to_file_with_motifs *path_to_HMDB_json_file* *path_to_store_spectrum_files*
*path_to_HMDB_mgf_file*
    path_to_file_with_motifs: a tab separated file with a selected motif, feature list and massql query on each line
    (output from make_pdf_with_smiles.py)
    path_to_HMDB_json_file: the path to the file with MS/MS spectra retrieved from GNPS in json format
    path_to_store_spectrum_files: folder where all the mgf-formatted text files with spectra will be stored that contain
    a selected motif (determined by MassQL)
    path_to_HMDB_mgf_file: the path to the file with mgf-formatted MS/MS spectra from GNPS
"""

# import statements
from massql import msql_engine
from sys import argv
from pathlib import Path
import re
import os
import pandas as pd
import sys, os

# functions
def parse_input(mgf_file: str) -> dict:
    """Parses mgf into strings, where each string is a spectrum and stores those in dict with the spectrum_id as key.

    :param mgf_file: str, name of mgf formatted file containing MS/MS spectra from GNPS
    :return: dictionary with {spectrum_id:record} where each record is a string containing the mgf-style accession of one
    compound
    """

    lines_mgf_file=open(mgf_file)
    spectrum_record_bool = False
    mgf_spectrum_records=[]
    mgf_spectrum_record = ""
    for line in lines_mgf_file:
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

    splitted_line = line.split("\t")
    assert len(
        splitted_line) == 3, "Expected a file with lines with tabs seperating motif, feature_list and massql query"
    motif, features, massql_query = splitted_line
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
        devnull = open(os.devnull, "w")

        # Suppress streams
        if self.suppress_stdout:
            self.original_stdout = sys.stdout
            sys.stdout = devnull

        if self.suppress_stderr:
            self.original_stderr = sys.stderr
            sys.stderr = devnull

    def __exit__(self, *args, **kwargs):
        # Restore streams
        if self.suppress_stdout:
            sys.stdout = self.original_stdout

        if self.suppress_stderr:
            sys.stderr = self.original_stderr

def search_motif_with_massql(massql_query: str, path_to_json_file: str) -> pd.DataFrame:
    """
    Uses MassQL to retrieve the identifiers of the spectra in the json file that contain the query

    :param massql_query: str, MassQL query made based on a Mass2Motif
    :param path_to_json_file: str, the path to a file with MS/MS spectra retrieved from GNPS in json format
    :return: returns a Pandas dataframe with the spectrum ids of the spectra as the index which contain the
    characteristics of the query.
    """
    # msql_engine only takes json files
    # the class Suppress suppresses the output to the command line that the msql_engine is giving through a print
    # function.
    with Suppress(suppress_stderr=True, suppress_stdout=True):
        df_massql_res=msql_engine.process_query(massql_query,path_to_json_file)
    if not df_massql_res.empty:
        df_massql_res.rename(columns={'scan': 'spectrum_id'}, inplace=True)
        df_massql_res.set_index("spectrum_id", inplace=True, drop=True)
        return df_massql_res

def make_mgf_file_for_spectra(motif: str, df_massql_matches: pd.DataFrame, path_to_store_spectrum_files: str,
                                       dict_with_mgf_spectra: dict) -> None:
    """
    Writes a spectrum containing the motif to a text file in mgf format using the identifier given by MassQL

    :param motif: str, the Mass2Motif for which the MassQL query was made and for which the spectrum was found, which
    contains the motif
    :param df_massql_matches: pd.Dataframe, with the spectrum ids of the spectra as the index which contain the
    characteristics of the query. The spectrum_id of GNPS spectra files formatted like this: CCMSLIB00000425029
    :param path_to_store_spectrum_files: str, folder where all the mgf-formatted text files (one for each motif) with
    spectra will be stored.
    :param dict_with_mgf_spectra: dict, with {spectrum_id:record} where each record is a string containing the mgf-style
    spectrum of one compound
    :return: path of the file containing the spectrum or spectra for one motif found by MassQL with the motif in the
    file name.
    """
    # the dataframe can be empty for some motifs, but if its not empty
    if df_massql_matches is not None:
        # create a path to store the mgf spectra of the matches found for the motif by MassQL
        path_to_spectra_file_for_motif = Path(
            fr"{path_to_store_spectrum_files}/mgf_spectra_for_{motif}_from_massql.txt")
        identifier_num = 0
        for identifier, row in df_massql_matches.iterrows():
            spectrum_record_mgf = dict_with_mgf_spectra[identifier]
            # if its the first identifier of the dataframe: make the mgf_spectrum file
            identifier_num += 1
            if identifier_num == 1:
                if os.path.exists(path_to_spectra_file_for_motif):
                    assert False, f"The mgf spectra file for {motif} already exists, remove it"
                # write to the new file
                else:
                    mgf_spectra_file = open(path_to_spectra_file_for_motif, "w")
                    mgf_spectra_file.write(spectrum_record_mgf)
                    mgf_spectra_file.write("\n")
                    mgf_spectra_file.close()
            # else append to the file with mgf spectra if the identifier is not the first identifier
            else:
                mgf_spectra_file = open(path_to_spectra_file_for_motif, "a")
                mgf_spectra_file.write("\n")
                mgf_spectra_file.write(spectrum_record_mgf)
                mgf_spectra_file.write("\n")
                mgf_spectra_file.close()
        return None

def main():
    """Main function of this module"""
    path_to_file_with_motifs = argv[1]
    path_to_HMDB_json_file = argv[2]
    path_to_store_spectrum_files = argv[3]
    path_to_mgf_file = argv[4]
    # step 1: parse the mgf file with spectra and put into dictionary with spectrum id's as keys.
    dict_with_mgf_spectra=parse_input(path_to_mgf_file)
    # step 2: parse the lines in the file where all the selected motifs and their corresponding massql queries are
    # listed.
    with open(path_to_file_with_motifs, "r") as lines_motif_file:
        for line in lines_motif_file:
            line = line.strip()
            line = line.replace('\n', '')
            # step 3: retrieve the motif and the massql query one by one
            motif, features, massql_query=parse_line_with_motif_and_query(line)
            # step 4: search for spectra with the motif in json file with MassQL
            df_massql_matches=search_motif_with_massql(massql_query, path_to_HMDB_json_file)
            # step 5: make a file with all the mgf formatted spectra for every identified match for a motif
            make_mgf_file_for_spectra(motif, df_massql_matches, path_to_store_spectrum_files,
                                                                       dict_with_mgf_spectra)

if __name__ == "__main__":
    main()