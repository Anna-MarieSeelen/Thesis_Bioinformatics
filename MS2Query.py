#!/usr/bin/env python3
"""
Author: Anna-Marie Seelen
Studentnumber:1008970
Description: alignment of protein sequences to a reference protein sequence
and returning stats
Usage: python3 *name_of_script* *path_folder_input* *path_folder_output*
    name_of_script: MS2Query.py
    path_folder_input: Define the folder in which your query spectra are stored. Accepted formats are: "mzML", "json",
    "mgf", "msp", "mzxml", "usi" or a pickled matchms object
    path_folder_output: Set the location where all your downloaded model files are stored
"""
# import statements
from sys import argv
import re
import subprocess
import os.path
from ms2query.run_ms2query import download_default_models, default_library_file_base_names, run_complete_folder
from ms2query.ms2library import create_library_object_from_one_dir

# functions
def ms2query(path_folder_input_files,path_folder_output_files):
    # Set the location where all your downloaded model files are stored
    ms2query_library_files_directory = "C:\\Users\\mseel\\OneDrive\\Documenten\\Dropbox\\Universiteit\\MBF derde jaar thesis\\MS2Query_output_test_files"
    # Define the folder in which your query spectra are stored.
    # Accepted formats are: "mzML", "json", "mgf", "msp", "mzxml", "usi" or a pickled matchms object.
    ms2_spectra_directory = "C:\\Users\\mseel\\OneDrive\\Documenten\\Dropbox\\Universiteit\\MBF derde jaar thesis\\test_file"

    # Downloads pretrained models and files for MS2Query (>10GB download)
    download_default_models(ms2query_library_files_directory, default_library_file_base_names())

    # Create a MS2Library object
    ms2library = create_library_object_from_one_dir(ms2query_library_files_directory, default_library_file_base_names())

    # Run library search and analog search on your files.
    run_complete_folder(ms2library, ms2_spectra_directory)

def main():
    """Main function of this module"""
    # step 1:
    path_folder_input_files = argv[1]
    path_folder_output_files = argv[2]

if __name__ == "__main__":
    main()