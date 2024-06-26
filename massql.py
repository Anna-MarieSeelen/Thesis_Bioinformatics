#!/usr/bin/env python3
"""
Author: Anna-Marie Seelen
Studentnumber:1008970
Description: Takes a file with MassQL queries, searches for spectra in json file that contain the queries and writes
these spectra to individual files in mgf-style.
Usage: python3 massql_search_spectra_with_motif.py *path_to_file_with_motifs_queries*
*pickle_file_with_gnps_all_positive_MS/MS_spectra* *path_to_store_spectrum_files* *path_to_store_match_files*
*path_to_store_mgf_file_and_name* *path_to_store_json_file_and_name*

    path_to_file_with_motifs_queries: a tab separated file with a selected motif, feature list and massql query on each
    line (output from make_pdf_with_smiles.py)
    pickle_file_with_gnps_all_positive_MS/MS_spectra: path to pickle file which is in:
    /mnt/LTR_userdata/hooft001/mass_spectral_embeddings/datasets/GNPS_15_12_21/ALL_GNPS_15_12_2021_positive_annotated.pickle
    path_to_store_spectrum_files: folder where all the mgf-formatted text files with spectra will be stored that contain
    a selected motif (determined by MassQL)
    path_to_store_match_files: folder where per Mass2Motif query a file will be stored containing the amount of massql
    matches and the smiles of each selected library spectrum
    path_to_store_mgf_file_and_name: the path and the file name that you want for the mgf-formatted MS/MS spectra from
    GNPS
    path_to_store_json_file_and_name: the path and the file name that you want for the json MS/MS spectra from
    GNPS
"""

# import statements
from massql import msql_engine
from sys import argv
from pathlib import Path
import re
import pandas as pd
import sys, os
from matchms import Scores, Spectrum
import json
from typing import List
import numpy
from pandas import json_normalize
import time
from matchms.filtering import select_by_relative_intensity
from matchms.exporting import save_as_mgf

# functions
def parse_input(mgf_file: str) -> dict:
    """Parses mgf into strings, where each string is a spectrum and stores those in dict with the spectrum_id as key.

    :param mgf_file: str, the path and the file name that you want for the mgf-formatted MS/MS spectra from
    GNPS
    :return: dictionary with {spectrum_id:record} where each record is a string containing the mgf-style accession of
    one compound
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

def make_json_mgf_file(path_to_pickle_file: str, path_to_store_json_file: str, path_to_store_mgf_file: str,
                           intensity_threshold=0.8) -> tuple:
    """
    Takes a pickle file with matchms spectrum objects and returns selected spectra in mgf and json format

    :param path_to_pickle_file: str, path to pickle file which contains matchms spectrum objects of GNPS library spectra
    :param path_to_store_json_file: str, the path and the file name that you want for the json MS/MS spectra from
    GNPS
    :param path_to_store_mgf_file: str, the path and the file name that you want for the mgf-formatted MS/MS spectra
    from GNPS
    :param intensity_threshold: float, relative intensity threshold to select fragment peaks to include in the json
    and mgf file. This threshold should be similar as the threshold used in the massql queries (default: 0.8).
    :return: tuple with strings, the strings are the path to the json and mgf file with the library spectra with the
    selected peaks
    """
    if os.path.exists(path_to_store_json_file):
        assert False, f"path to json file {path_to_store_json_file} exists, remove it!"
    if os.path.exists(path_to_store_mgf_file):
        assert False, f"path to mgf file {path_to_store_mgf_file} exists, remove it!"
    obj = pd.read_pickle(path_to_pickle_file)
    spectrum_list=[]
    for spectrum in obj:
        spectrum = select_by_relative_intensity(spectrum, intensity_from=intensity_threshold)
        spectrum_list.append(spectrum)
    #there is some weird error that save_as_json gives but ignore it
    save_as_json(spectrum_list,path_to_store_json_file)
    save_as_mgf(spectrum_list, path_to_store_mgf_file)
    return path_to_store_json_file, path_to_store_mgf_file

def save_as_json(spectrums: List[Spectrum], filename: str):
    """Save spectrum(s) as json file.
    :py:attr:`~matchms.Spectrum.losses` of spectrum will not be saved.
    Arguments:
    ----------
    spectrums:
        Expected input are match.Spectrum.Spectrum() objects.
    filename:
        Provide filename to save spectrum(s).

    Note: copied and adjusted function from matchms
    """
    if not isinstance(spectrums, list):
        # Assume that input was single Spectrum
        spectrums = [spectrums]

    # Write to json file
    with open(filename, 'w') as fout:
        fout.write("[")
        for i, spectrum in enumerate(spectrums):
            spec = spectrum.clone()
            peaks_list = str(numpy.vstack((spec.peaks.mz, spec.peaks.intensities)).T.tolist())

            # Convert matchms.Spectrum() into dictionaries
            spectrum_dict = {key: spec.metadata[key] for key in spec.metadata}
            spectrum_dict["spectrum_id"] = spectrum_dict["spectrumid"]
            del spectrum_dict["spectrumid"]
            spectrum_dict["peaks_json"] = peaks_list
            spectrum_dict["Precursor_MZ"] = spectrum_dict["precursor_mz"]
            del spectrum_dict["precursor_mz"]
            spectrum_dict["scan"] = spectrum_dict["scans"]
            del spectrum_dict["scans"]

            json.dump(spectrum_dict, fout)
            if i < len(spectrums) - 1:
                fout.write(",")
        fout.write("]")

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

def read_json(json_file: str) -> pd.DataFrame:
    """
    Loads a json file into a pd.DataFrame

    :param json_file: str, path to json file with selected GNPS library spectra in json format
    :return: pd.DataFrame with GNPS library spectra with selected peaks in json format
    """
    with open(json_file, 'r') as f:
        dict=json.load(f)
    df_json=json_normalize(dict)
    df_json.set_index("spectrum_id",inplace=True, drop=True)
    return df_json

def search_motif_with_massql(massql_query: str, path_to_json_file: str) -> pd.DataFrame:
    """
    Uses MassQL to retrieve the identifiers of the spectra in the json file that contain the query

    :param massql_query: str, MassQL query made based on a Mass2Motif
    :param path_to_json_file: str, the path to a file with MS/MS spectra retrieved from GNPS in json format
    :return: returns a Pandas dataframe with the spectrum ids of the spectra as the index which contain the
    characteristics of the query.
    """
    # msql_engine only takes json files and will print search bars, just ignore it.
    df_massql_res=msql_engine.process_query(massql_query,path_to_json_file)
    if not df_massql_res.empty:
        df_massql_res.rename(columns={'scan': 'spectrum_id'}, inplace=True)
        df_massql_res.set_index("spectrum_id", inplace=True, drop=True)
        return df_massql_res

def select_massql_matches(path_to_store_match_files: str, df_massql_matches: pd.DataFrame, df_json: pd.DataFrame,
                              motif: str) -> pd.DataFrame:
    """
    Makes a dataframe and csv file with 30 or less matched library spectra for a motif query and the precmz, and smiles

    :param path_to_store_match_files: folder where per Mass2Motif query a file will be stored containing the amount of
    massql matches and the smiles of each selected library spectrum
    :param df_massql_matches: pd.DataFrame, with the spectrum ids of the spectra as the index which contain the
    characteristics of the query.
    :param df_json: pd.DataFrame, with GNPS library spectra with selected peaks in json format
    :param motif: str, the mass2motif for which the query was made
    :return: pd.DataFrame with the selected spectrum ids as index, the precmz of the selected spectra and the smiles
    """

    df_matched_spectra_smiles = pd.merge(df_massql_matches["precmz"], df_json[["Precursor_MZ", "smiles"]],
                                         left_index=True, right_index=True)
    # for some matches there are no smiles so remove those from the dataframe
    df_matched_spectra_smiles.drop(df_matched_spectra_smiles.index[df_matched_spectra_smiles['smiles'] == 'N/A'],
                                   inplace=True)
    df_matched_spectra_smiles.drop(df_matched_spectra_smiles.index[df_matched_spectra_smiles['smiles'] == ' '],
                                   inplace=True)
    df_matched_spectra_smiles.drop_duplicates('smiles', inplace=True)
    df_matched_spectra_smiles.reset_index(inplace=True)
    df_matched_spectra_smiles=df_matched_spectra_smiles.set_index('spectrum_id')
    # if you have more than 30 non-redundant matches randomly select 30 matches, because otherwise it will take too long
    if len(df_matched_spectra_smiles) > 30:
        df_matched_spectra_smiles=df_matched_spectra_smiles.sample(n=30)
    count=len(df_matched_spectra_smiles)
    df_matched_spectra_smiles.to_csv(f"{path_to_store_match_files}/{motif}_amount_of_matches_{count}", sep='\t',
                                     encoding='utf-8')
    return df_matched_spectra_smiles

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
    :return: None
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
    before_script = time.perf_counter()
    path_to_file_with_motifs_queries = argv[1]
    path_to_pickle_file = argv[2]
    path_to_store_spectrum_files = argv[3]
    path_to_store_match_files=argv[4]
    path_to_store_mgf_file = argv[5]
    path_to_store_json_file=argv[6]
    # step 1: transfer the data from the pickle file to a json and mgf file and select spectrum peaks above threshold
    path_to_store_json_file, path_to_store_mgf_file = make_json_mgf_file(path_to_pickle_file, path_to_store_json_file,
                                                                         path_to_store_mgf_file,
                                                                         intensity_threshold=0.8)
    # step 2: make a dictionary with all the accessions in the mgf file
    dict_with_mgf_spectra=parse_input(path_to_store_mgf_file)
    # step 3: read the json file into a dataframe
    df_json=read_json(path_to_store_json_file)
    # step 4: parse the lines in the file where all the selected motifs and their corresponding massql queries are
    # listed.
    with open(path_to_file_with_motifs_queries, "r") as lines_motif_file:
        for line in lines_motif_file:
            line = line.strip()
            line = line.replace('\n', '')
            # step 5: retrieve the motif and the massql query one by one
            motif, features, massql_query=parse_line_with_motif_and_query(line)
            # step 6: search for spectra with the motif in json file with MassQL
            df_massql_matches=search_motif_with_massql(massql_query, path_to_store_json_file)
            # step 7: make a df and csv file with the smiles of the selected library spectra for each motif
            df_massql_matches_with_smiles = select_massql_matches(path_to_store_match_files, df_massql_matches, df_json,
                                                                  motif)
            # step 8: make a file with all the mgf formatted library spectra for every identified match for a motif
            make_mgf_file_for_spectra(motif, df_massql_matches_with_smiles, path_to_store_spectrum_files,
                                                                       dict_with_mgf_spectra)
    after_script = time.perf_counter()
    print("how long the total script took {0}".format(after_script - before_script))

if __name__ == "__main__":
    main()
