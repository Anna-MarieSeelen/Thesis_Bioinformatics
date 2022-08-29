#!/usr/bin/env python3
"""
Author: Anna-Marie Seelen
Studentnumber:1008970
Description: run MS2Query
Usage: python3 *name_of_script* *path_folder_input* *path_folder_library_files*
    name_of_script: MS2Query.py
    path_folder_input: Define the folder in which your query spectra are stored. Accepted formats are: "mzML", "json",
    "mgf", "msp", "mzxml", "usi" or a pickled matchms object. The results generated by MS2Query, are stored as csv files
    in a results directory within the same directory as your query spectra.
    path_folder_library_files: Set the location where all your downloaded MS2Query libraries are stored
"""
# import statements
from sys import argv
import re
import subprocess
import os.path
from ms2query.run_ms2query import download_default_models, default_library_file_base_names, run_complete_folder
from ms2query.ms2library import create_library_object_from_one_dir

# functions
def ms2query(path_folder_input_files: str,path_folder_library_files: str) -> None:
    """
    Downloads GNPS mass spectral library, executes MS2Query, and makes folder called where output csv file is stored.
    :param path_folder_input_files: str, path of location of folder in which your query spectra are stored
    :param path_folder_library_files: str, path of location of folder where all your downloaded MS2Query libraries are
    stored
    :return: no return, but makes folder with csv file containing matches and analogues for every compound in input_file
    """

    # Downloads pretrained models and files for MS2Query (>10GB download)
    download_default_models(path_folder_library_files, default_library_file_base_names())
    # Create a MS2Library object
    ms2library = create_library_object_from_one_dir(path_folder_library_files, default_library_file_base_names())
    # Run library search and analog search on your files.
    run_complete_folder(ms2library, path_folder_input_files)

def main() -> None:
    """Main function of this module"""
    # step 1: execute MS2Query
    path_folder_input_files = argv[1]
    path_folder_library_files = argv[2]
    ms2query(path_folder_input_files, path_folder_library_files)

if __name__ == "__main__":
    main()