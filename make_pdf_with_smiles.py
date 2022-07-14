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
def connection_Mass2Motifs_to_documents(path_file_with_MS2LDA_csv):

    #read both csv files into dataframe
    df_MS2LDA = pd.read_csv(path_file_with_MS2LDA_csv, header=0)

    #add Motif and Document number into 1 dataframe, but only include Motifs with more than 5 doc --> you will lose some doc
    df_MS2LDA["Document+Probability+Overlap"] = df_MS2LDA[["Document", "Probability", "Overlap Score"]].values.tolist()
    df_new2 = df_MS2LDA[["Motif", "Document+Probability+Overlap"]].groupby("Motif", as_index=True).aggregate({"Document+Probability+Overlap":list})
    df_new2 = df_new2[df_new2["Document+Probability+Overlap"].str.len() >= 5]
    df_new=df_MS2LDA[["Motif","Document"]].groupby("Motif", as_index=True).aggregate({"Document":list})
    # we only want the motifs found in more than 5 spectra
    df_new=df_new[df_new['Document'].str.len() >= 5]
    # solution from: https://stackoverflow.com/questions/58297277/how-to-get-a-length-of-lists-in-pandas-dataframe
    f = lambda x: 'document_{}'.format(x + 1)
    df= pd.DataFrame(df_new.Document.values.tolist(),df_new.index, dtype=object).fillna('').rename(columns=f)
    # solution from: https://stackoverflow.com/questions/44663903/pandas-split-column-of-lists-of-unequal-length-into-multiple-columns
    df = pd.concat([df, df_new2], axis=1)
    return df

def information_document_node(path_to_file_with_MS2Query_csv):
    df_MS2Query = pd.read_csv(path_to_file_with_MS2Query_csv, header=0)
    df_sub = df_MS2Query[df_MS2Query.ms2query_model_prediction > 0.7]
    df_part=pd.DataFrame(df_sub[["query_spectrum_nr", "precursor_mz_difference", "smiles" ]]).set_index("query_spectrum_nr")
    df_new = pd.DataFrame(df_MS2Query[["query_spectrum_nr","precursor_mz_query_spectrum", "ms2query_model_prediction", "precursor_mz_analog"]]).set_index("query_spectrum_nr")
    df=pd.concat([df_new, df_part], axis=1).fillna('')
    return df

#Make the other file with fragements input too: Mass2Motif, fragments
def information_document_Mass2Motif_node(path_file_with_MS2LDA_csv_fragments):
    df_Motif_fragments= pd.read_csv(path_file_with_MS2LDA_csv_fragments, header=0)
    df_Motif_fragments["Fragment+Probability"]=df_Motif_fragments[["Feature", "Probability"]].values.tolist()
    #.sort(key = lambda x: x[1], reverse=True)
    #https: // www.geeksforgeeks.org / python - sort - list - according - second - element - sublist /
    df_new = df_Motif_fragments[["Motif", "Fragment+Probability"]].groupby("Motif", as_index=True).aggregate({"Fragment+Probability": list})
    # TODO: adjust sorted thing here
    #list.sort(key=lambda x: x[1], reverse=True
    # to defide the list of lists into more columns
    df_new2 = df_new["Fragment+Probability"]
    f = lambda x: 'fragment_{}+probability'.format(x + 1)
    df_new=df_new["Fragment+Probability"].apply(pd.Series).fillna('').rename(columns=f)
    # https://stackoverflow.com/questions/45107523/pandas-convert-list-of-lists-to-multiple-columns
    df_new=pd.concat([df_new, df_new2], axis=1)
    return df_new

def visualize_mol(smiles):
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

def main() -> None:
    """Main function of this module"""
    path_to_file_with_MS2Query_csv= argv[1]
    path_file_with_MS2LDA_csv = argv[2]
    path_file_with_Motif_fragments_csv = argv[3]
    # step 1: rearrange MS2LDA file with mass2motifs and documents in dataframe
    df_motifs_to_doc=connection_Mass2Motifs_to_documents(path_file_with_MS2LDA_csv)
    #step 2: rearrange MS2Query file with documents and smiles of analogs in dataframe
    df_doc_to_smiles=information_document_node(path_to_file_with_MS2Query_csv)
    # step 3: rearrange MS2LDA fragment file with mass2motifs and fragments in dataframe
    df_motifs_to_frag=information_document_Mass2Motif_node(path_file_with_Motif_fragments_csv)
    #step 4: make a PDF where you can see each Mass2Motif, its fragments, and the associated analog structures
    pdf=FPDF()
    pdf.add_page()
    pdf.set_font("helvetica", size=10)
    for index, row in df_motifs_to_doc.iterrows():
        amount_of_smiles = 0
        # calculating the amount of smiles with a probability > 0.7 associated with each motif
        for cell in row:
            if cell in df_doc_to_smiles.index.values.tolist():
                if df_doc_to_smiles.at[cell, "smiles"] != "":
                    amount_of_smiles+=1
        # selection of printed Mass2Motifs based on the amount of smiles with a probability > 0.7. At least 2 or more Smiles are needed.
        if amount_of_smiles >= 2:
            pdf.multi_cell(200,5, txt="\n", align = 'C')
            pdf.set_font("helvetica", "B", size=10)
            # writing the mass_2_motif name, which is the index of df_motifs_to_doc in the file
            pdf.multi_cell(200, 5, txt="{0}\n".format(index), align = 'C')
            pdf.set_font("helvetica", size=10)
            # if the mass2motif is in the list of indexes of df_motifs_to_frag then look for the column with fragments + probability of that mass2motif and sort that list on the probability value
            if index in df_motifs_to_frag.index.values.tolist():
                pdf.multi_cell(200, 5, txt="{0}\n".format(sorted(df_motifs_to_frag.at[index, "Fragment+Probability"], key = lambda x: x[1], reverse=True)), align = 'L')
            else:
                pdf.multi_cell(200, 5, txt="motif does not have fragments with probability > 0.05", align = 'L')
            # print the "Document+probability+Overlap" column for the selected index, which is the selected Mass2Motif
            pdf.set_font("helvetica", "I", size=10)
            pdf.multi_cell(200, 5, txt="[[Spectrum number, probability, overlap score], [Spectrum number, probability, overlap score], etc.] :\n", align='L')
            pdf.set_font("helvetica", size=10)
            pdf.multi_cell(200, 5, txt="{0}\n".format(df_motifs_to_doc.at[index, "Document+Probability+Overlap"]), align = 'L')
            # print the amount of lists in the "Document+probability+Overlap" column for the selected index, which is equal to the amount of documents the motif is associated with
            pdf.multi_cell(200, 5, txt="amount of spectra associated with motif: {0}\n".format(len(df_motifs_to_doc.at[index, "Document+Probability+Overlap"])), align = 'L')
            # for each individual cell in df_motifs_to_doc, which is the document number associated with the index, so the mass2motif
            for cell in row:
                # if the cell value is an int, so if its really a spectrum number (because there is also a list of lists in this dataframe)
                if isinstance(cell, int):
                    # if the spectrum number is in the index of df_doc_to_smiles
                    if cell in df_doc_to_smiles.index.values.tolist():
                        # if the smiles associated with the spectrum number is not empty, because most cells are empty because of a too low prediction probability of the analogue
                        if df_doc_to_smiles.at[cell, "smiles"] != "":
                            # then print the associated spectrum number that the smiles belongs to and the precursor m/z of the spectrum
                            pdf.multi_cell(200, 10, txt="spectrum number: {0}, precursor m/z: {1}\n".format(cell, df_doc_to_smiles.at[cell, "precursor_mz_query_spectrum"]),align='L')
                            # also use the fuction visualize_mol to print a picture of the structure of the analog using the smiles that is associated with the spectrum
                            pdf.image(visualize_mol(df_doc_to_smiles.at[cell, "smiles"]))
                            # print the predicition probability of the analog, and the mz of the smiles that is associated with the spectrum
                            pdf.multi_cell(200, 10, txt="analog m/z:{0}, prediction: {1}\n".format(
                                df_doc_to_smiles.at[cell, "precursor_mz_analog"],
                                df_doc_to_smiles.at[cell, "ms2query_model_prediction"]), align='L')
    pdf.output("output.pdf")

if __name__ == "__main__":
    main()