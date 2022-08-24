#!/usr/bin/env python3
"""
Author: Anna-Marie Seelen
Studentnumber:1008970
Description: run MS2Query
Usage: python3 *name_of_script* *path_to_file_with_MS2Query_csv* *path_file_with_MS2LDA_csv*
    name_of_script: make_pdf_with_smiles.py
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
import rdkit.Chem as Chem
from rdkit.Chem.Draw import MolToImage
from fpdf import FPDF

# functions
def connection_Mass2Motifs_to_documents(path_to_file_with_MS2Query_csv, path_file_with_MS2LDA_csv):

    #read both csv files into dataframe
    df_MS2Query = pd.read_csv(path_to_file_with_MS2Query_csv, header=0)
    df_MS2LDA = pd.read_csv(path_file_with_MS2LDA_csv, header=0)

    #add Motif and Document number into 1 dataframe, but only include Motifs with more than 5 doc --> you will lose some doc
    df_new=df_MS2LDA[["Motif","Document"]].groupby("Motif", as_index=True).aggregate({"Document":list})
    # we only want the motifs found in more than 5 spectra
    df_new=df_new[df_new['Document'].str.len() >= 5]
    df_new2 = df_new["Document"]
    #https://stackoverflow.com/questions/58297277/how-to-get-a-length-of-lists-in-pandas-dataframe
    f = lambda x: 'document_{}'.format(x + 1)
    df= pd.DataFrame(df_new.Document.values.tolist(),df_new.index, dtype=object).fillna('').rename(columns=f)
    # https://stackoverflow.com/questions/44663903/pandas-split-column-of-lists-of-unequal-length-into-multiple-columns
    df = pd.concat([df, df_new2], axis=1)
    return df

def information_document_node(path_to_file_with_MS2Query_csv):
    df_MS2Query = pd.read_csv(path_to_file_with_MS2Query_csv, header=0)
    df_sub = df_MS2Query[df_MS2Query.ms2query_model_prediction > 0.7]
    df_part=pd.DataFrame(df_sub[["query_spectrum_nr","ms2query_model_prediction", "precursor_mz_difference", "smiles" ]]).set_index("query_spectrum_nr")
    df_new = pd.DataFrame(df_MS2Query[["query_spectrum_nr","precursor_mz_query_spectrum"]]).set_index("query_spectrum_nr")
    df=pd.concat([df_new, df_part], axis=1).fillna('')
    return df

#Make the other file with fragements input too: Mass2Motif, fragments
def information_document_Mass2Motif_node(path_file_with_MS2LDA_csv_fragments):
    df_Motif_fragments= pd.read_csv(path_file_with_MS2LDA_csv_fragments, header=0)
    df_Motif_fragments["Fragment+Probability"]=df_Motif_fragments[["Feature", "Probability"]].values.tolist()
    print(df_Motif_fragments)
    #.sort(key = lambda x: x[1], reverse=True)
    #https: // www.geeksforgeeks.org / python - sort - list - according - second - element - sublist /
    df_new = df_Motif_fragments[["Motif", "Fragment+Probability"]].groupby("Motif", as_index=True).aggregate({"Fragment+Probability": list})
    print(df_new)
    # to defide the list of lists into more columns
    df_new2 = df_new["Fragment+Probability"]
    print(df_new2)
    f = lambda x: 'fragment_{}+probability'.format(x + 1)
    df_new=df_new["Fragment+Probability"].apply(pd.Series).fillna('').rename(columns=f)
    # https://stackoverflow.com/questions/45107523/pandas-convert-list-of-lists-to-multiple-columns
    df_new=pd.concat([df_new, df_new2], axis=1)
    return df_new

def make_df_smiles(df_doc_to_smiles):
    #TODO: waarom doet deze smiles het niet?
    data_smiles=[]
    for index,row in df_doc_to_smiles.iterrows():
        if row["smiles"] != "":
            data_smiles.append(row)
        else:
            pass
    df_smiles = pd.DataFrame(data_smiles)
    return df_smiles

def visualize_mol(smiles):
    #TODO: doesn't work for some reason: cannot put it in pdf
    """
    function adapted from: https://pchanda.github.io/See-substructure-in-molecule/
    :param smiles:
    :return:
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    Chem.Kekulize(mol)
    img = MolToImage(mol, size=(200, 200), fitImage=True)
    return img

def visualize_on_command_line(img):
    """Run jellyfish program on fasta file
    input_fn: string, filename of input FASTA file
    kmer_size: int, size of k-mers used by jellyfish
    """
    #out_fn = 'tomato{}'.format(kmer_size)
    cmd = 'viu {}' \
            .format(img)
    # e = subprocess.check_output(cmd, shell=True)
    # if os.path.exists(out_fn):
    #     cmd = 'jellyfish stats {}'.format(out_fn)
    # else:
    #     cmd = 'jellyfish stats {}_0'.format(out_fn)
    res = subprocess.check_output(cmd, shell=True)
    return res

def main() -> None:
    """Main function of this module"""
    # step 1: parse MS2Query file
    path_to_file_with_MS2Query_csv= argv[1]
    # step 2: parse MS2LDA file
    path_file_with_MS2LDA_csv = argv[2]
    # step 3: input Mass2Motif fragments file
    path_file_with_Motif_fragments_csv = argv[3]
    # stept 3: put relevant info in new document
    df_motifs_to_doc=connection_Mass2Motifs_to_documents(path_to_file_with_MS2Query_csv, path_file_with_MS2LDA_csv)
    df_doc_to_smiles=information_document_node(path_to_file_with_MS2Query_csv)
    df_motifs_to_frag=information_document_Mass2Motif_node(path_file_with_Motif_fragments_csv)
    # step 5: make dataframe with only smiles
    df_smiles=make_df_smiles(df_doc_to_smiles)

if __name__ == "__main__":
    main()