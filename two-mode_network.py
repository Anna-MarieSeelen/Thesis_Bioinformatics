#!/usr/bin/env python3
"""
Author: Anna-Marie Seelen
Studentnumber:1008970
Description: run MS2Query
Usage: python3 *name_of_script* *path_to_file_with_MS2Query_csv* *path_file_with_MS2LDA_csv*
    name_of_script: two-mode_network.py
    path_to_file_with_MS2Query_csv: Define the folder in which your query spectra are stored. Accepted formats are: "mzML", "json",
    "mgf", "msp", "mzxml", "usi" or a pickled matchms object. The results generated by MS2Query, are stored as csv files
    in a results directory within the same directory as your query spectra.
    path_file_with_MS2LDA_csv: Set the location where all your downloaded MS2Query libraries are stored
"""
# import statements
from sys import argv
import pandas as pd
import re
import subprocess
import os.path

# functions
def connection_Mass2Motifs_to_documents(path_to_file_with_MS2Query_csv, path_file_with_MS2LDA_csv):

    #read both csv files into dataframe
    df_MS2Query = pd.read_csv(path_to_file_with_MS2Query_csv, header=0)
    df_MS2LDA = pd.read_csv(path_file_with_MS2LDA_csv, header=0)

    #add Motif and Document number into 1 dataframe
    df_new=df_MS2LDA[["Motif","Document"]].groupby("Motif", as_index=True).aggregate({"Document":list})
    f = lambda x: 'document_{}'.format(x + 1)
    df= pd.DataFrame(df_new.Document.values.tolist(),df_new.index, dtype=object).fillna('').rename(columns=f)
    #https://stackoverflow.com/questions/44663903/pandas-split-column-of-lists-of-unequal-length-into-multiple-columns
    return df

def information_document_node(path_to_file_with_MS2Query_csv):
    df_MS2Query = pd.read_csv(path_to_file_with_MS2Query_csv, header=0)
    df_sub = df_MS2Query[df_MS2Query.ms2query_model_prediction > 0.7]
    df_part=pd.DataFrame(df_sub[["query_spectrum_nr","ms2query_model_prediction", "precursor_mz_difference", "smiles" ]]).set_index("query_spectrum_nr")
    df_new = pd.DataFrame(df_MS2Query[["query_spectrum_nr","precursor_mz_query_spectrum"]]).set_index("query_spectrum_nr")
    df=pd.concat([df_new, df_part], axis=1).fillna('')
    return df

#Make the other file with fragements input too: Mass2Motif, fragments
def information_document_Mass2Motif_node(path_file_with_MS2LDA_csv):
    return None

def to_output_tab_file(alignment_list): #this should be csv format, not done yet
    """Takes a nested sorted list and outputs a tab delimited file

    alignment_list: nested list with families and alignment lenghts
    return: tab delimited text file with the contents of each sub list on a line
    """
    tab_file=open("out_tab_delimited.txt", "w")
    for sub_list in alignment_list:
        tab_file.write("{0}\t{1}\t{2}".format(sub_list[0], sub_list[1], sub_list[2]))
        tab_file.write("\n")
    return None

def main() -> None:
    """Main function of this module"""
    # step 1: parse MS2Query file
    path_to_file_with_MS2Query_csv= argv[1]
    # step 2: parse MS2LDA file
    path_file_with_MS2LDA_csv = argv[2]
    # stept 3: put relevant info in new document
    #print(connection_Mass2Motifs_to_documents(path_to_file_with_MS2Query_csv, path_file_with_MS2LDA_csv))
    print(information_document_node(path_to_file_with_MS2Query_csv))

if __name__ == "__main__":
    main()