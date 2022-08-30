#!/usr/bin/env python3
"""
Author: Anna-Marie Seelen
Studentnumber:1008970
Description: takes result csv files of MS2LDA and MS2Query and outputs a pdf that gives an overview of
motifs and analogs. In addition, it outputs a text file with the selected mass2motifs and their MassQL querries.
Usage: python3 *name_of_script* *path_to_file_with_MS2Query_csv* *path_file_with_MS2LDA_csv*
*path_file_with_Motif_fragments_csv*

    name_of_script: make_pdf_with_smiles.py
    path_to_file_with_MS2Query_csv: output file from the MS2Query.py script with document number and
    information about the matched compound
    path_file_with_MS2LDA_csv: file downloaded from MS2LDA.org which contains document number and name of the motif
    associated with the document
    path_file_with_Motif_fragments_csv: file downloaded from MS2LDA.org with contains the fragments and neutral losses
    and which motif these are associated with
"""

# import statements
from sys import argv
import pandas as pd
import rdkit.Chem as Chem
from rdkit.Chem.Draw import MolToImage
from fpdf import FPDF
import re
from pathlib import Path
import os

# functions
def convert_Mass2Motifs_to_df(path_file_with_MS2LDA_csv: str, threshold_spectra_in_motif=5) -> pd.DataFrame:
    """
    Takes a csv file from MS2LDA.org and rearranges the motifs and documents and returns them in a pandas dataframe.

    :param path_file_with_MS2LDA_csv: str, path to location of file with csv from MS2LDA containing motifs and documents
    :param threshold_spectra_in_motif: int, default=5, minimum amount of spectra associated with a motif for the motif
     to be included in the output dataframe
    :return: pandas dataframe with the motifs as index and each associated document in a separate column. The last
    column is a list of lists of the associated document, the probability and the overlap score.
    """

    df_MS2LDA = pd.read_csv(path_file_with_MS2LDA_csv, header=0)
    #add a new column the dataframe with a list of the document, probability and overlap score for each document
    df_MS2LDA["Document+Probability+Overlap"] = df_MS2LDA[["Document", "Probability", "Overlap Score"]].values.tolist()
    # First make a dataframe with the Motifs as index and a second column with a list of lists containing the associated
    # documents, probability and overlap scores
    df_motif_and_list = df_MS2LDA[["Motif", "Document+Probability+Overlap"]].groupby("Motif", as_index=True).aggregate(
        {"Document+Probability+Overlap": list})
    # We only want the motifs found in more than x spectra in the dataframe
    df_motif_and_list = df_motif_and_list[df_motif_and_list["Document+Probability+Overlap"].str.len() >= threshold_spectra_in_motif]
    # Make a second dataframe with the Motifs as index and a second column with a list containing the associated
    # documents
    df_motif_and_doclist = df_MS2LDA[["Motif","Document"]].groupby("Motif", as_index=True).aggregate({"Document":list})
    # Again we only want the motifs found in more than x spectra in the dataframe
    df_motif_and_doclist = df_motif_and_doclist[df_motif_and_doclist['Document'].str.len() >= threshold_spectra_in_motif]
    # solution from: https://stackoverflow.com/questions/58297277/how-to-get-a-length-of-lists-in-pandas-dataframe
    # now we want each document to be in a separate column in the dataframe
    f = lambda x: 'document_{}'.format(x + 1)
    df_motif_and_doc = pd.DataFrame(df_motif_and_doclist.Document.values.tolist(), df_motif_and_doclist.index,
                                    dtype=object).fillna('').rename(columns=f)
    # solution from:
    # https://stackoverflow.com/questions/44663903/pandas-split-column-of-lists-of-unequal-length-into-multiple-columns
    df_motifs_and_doc_and_list = pd.concat([df_motif_and_doc, df_motif_and_list], axis=1)
    return df_motifs_and_doc_and_list

def convert_ms2query_to_df(path_to_file_with_MS2Query_csv: str, threshold_ms2query_pred=0.7) -> pd.DataFrame:
    """
    Takes resulting csv file from MS2Query, selects the documents and analogs, returns them in a pandas dataframe.

    :param path_to_file_with_MS2Query_csv: str, path to location of file with result csv from MS2Query containing
    documents and their analogs
    :param threshold_ms2query_pred: int, default=0.7, minimal ms2query prediction for the analogue to be included in the
    output dataframe
    :return: pandas dataframe with the associated document and their precursor m/z. Prediction certainty,
    precursor_mz_analog, and smiles of analogs are given for prediction certainties above threshold.
    """

    df_MS2Query = pd.read_csv(path_to_file_with_MS2Query_csv, header=0)
    #selecting the spectra numbers where the ms2query prediction is above the threshold
    df_analog_above_threshold = df_MS2Query[df_MS2Query.ms2query_model_prediction > threshold_ms2query_pred]
    #making a dataframe where the selected the columns can go into
    df_analog_above_threshold_and_info = pd.DataFrame(df_analog_above_threshold[
                                                          ["query_spectrum_nr", "ms2query_model_prediction",
                                                           "precursor_mz_analog", "precursor_mz_difference",
                                                           "smiles"]]).set_index("query_spectrum_nr")
    #making the second dataframe with all spectra numbers and corresponding precursor masses
    df_all_spectra_and_mz = pd.DataFrame(df_MS2Query[["query_spectrum_nr", "precursor_mz_query_spectrum"]]).set_index(
        "query_spectrum_nr")
    df_all_spectra_and_smiles_above_threshold = pd.concat([df_all_spectra_and_mz, df_analog_above_threshold_and_info],
                                                          axis=1).fillna('')
    return df_all_spectra_and_smiles_above_threshold

def convert_fragments_in_motif_to_df(path_file_with_MS2LDA_csv_fragments: str) -> pd.DataFrame:
    """
    Takes a csv file from MS2LDA.org and rearranges the motifs and features and returns them in a pandas dataframe.

    :param path_file_with_MS2LDA_csv_fragments: str, path to location of csv file from MS2LDA containing motifs and
    features
    :return: pandas dataframe with the motifs as index and each feature and probability in a list in separate column.
    The last column contains a lists of list with each associated fragment and probability
    """

    df_Motif_fragments= pd.read_csv(path_file_with_MS2LDA_csv_fragments, header=0)
    # First make a column with a list of the probability and feature combined
    df_Motif_fragments["Fragment+Probability"]=df_Motif_fragments[["Feature", "Probability"]].values.tolist()
    # First make a dataframe with the Motifs as index and a second column with a list of lists containing the fragments,
    # and probability
    df_motif_list_of_lists_feature = df_Motif_fragments[["Motif", "Fragment+Probability"]].groupby("Motif", as_index=True).aggregate({"Fragment+Probability": list})
    # Make a second dataframe with the Motifs as index and a second column with a list containing the associated
    # features and probability in a separate column
    f = lambda x: 'fragment_{}+probability'.format(x + 1)
    df_motif_and_sep_fragments = df_motif_list_of_lists_feature["Fragment+Probability"].apply(pd.Series).fillna(
        '').rename(columns=f)
    # https://stackoverflow.com/questions/45107523/pandas-convert-list-of-lists-to-multiple-columns
    df_new=pd.concat([df_motif_and_sep_fragments, df_motif_list_of_lists_feature], axis=1)
    return df_new

def visualize_mol(smiles: str):
    """
    Takes a smiles as input and outputs an PIL<PNG> image

    :param smiles: str, a smiles string that corresponds to a molecular structure
    :return: PIL<PNG> image, image of structure corresponding with the smiles

    function adapted from: https://pchanda.github.io/See-substructure-in-molecule/
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    Chem.Kekulize(mol)
    img = MolToImage(mol, size=(200, 200), fitImage=True)
    return img

def make_MassQL_search(fragments: list) -> str:
    """
    Takes a list of lists containing fragments and probabilities and returns corresponding MassQL query

    :param fragments: list of lists containing fragments and probabilities belonging to a Mass2Motif
    :return: str, the corresponding MassQL query
    """
    query="QUERY scaninfo(MS2DATA) WHERE POLARITY = Positive " #MS1DATA doesn't give any results, so MS2DATA it is
    fragments=sorted(fragments, key = lambda x: x[1], reverse=True)
    print(fragments)
    fragments_rel_intensity=list(map(list, fragments)) # to make a copy, but not a reference of fragments
    for i in range(len(fragments_rel_intensity)):
        if i==0:
            fragments_rel_intensity[i][1]=1.0
        else:
            fragments_rel_intensity[i][1]=round(fragments[i][1]/fragments[0][1],3)
    print(fragments_rel_intensity)
    for fragment in fragments_rel_intensity:
        for string in fragment:
            if fragment[1]==1.0:
                if type(string) == str:
                    if re.search(r'loss', string) != None:
                        query += "AND MS2NL = {0}:TOLERANCEMZ={1}:INTENSITYMATCH=Y:INTENSITYMATCHREFERENCE ".format(re.search(r'\_(.*)', string).group(1), 0.01)
                    else:
                        query += "AND MS2PROD = {0}:TOLERANCEMZ={1}:INTENSITYMATCH=Y:INTENSITYMATCHREFERENCE ".format(re.search(r'\_(.*)', string).group(1),
                                                                             0.01)
                else:
                    pass
            else:
                if type(string)==str:
                    if re.search(r'loss', string) != None:
                        query+="AND MS2NL = {0}:TOLERANCEMZ={1}:INTENSITYMATCH=Y*{2}:INTENSITYMATCHPERCENT={3} ".format(re.search(r'\_(.*)', string).group(1), 0.01,fragment[1],20)
                    else:
                        query+="AND MS2PROD = {0}:TOLERANCEMZ={1}:INTENSITYMATCH=Y*{2}:INTENSITYMATCHPERCENT={3} ".format(re.search(r'\_(.*)', string).group(1), 0.01,fragment[1],20)
                else:
                    pass
    print(query)
    return query

def make_file_with_massql_querries(df_motifs_to_frag: pd.DataFrame, list_of_selected_motifs: list) -> str:
    """Outputs a tab delimited file with the selected motifs, their fragments and their MassQL querries

    :param df_motifs_to_frag: pandas dataframe, motifs as index and each feature and probability in a list in
    separate column. The last column contains a lists of list with each associated fragment and probability
    :param list_of_selected_motifs: list, motifs found in more than 5 spectra in the dataframe and that contain at least
    one fragments/neutral losses with a prob > 0.05 and at least 2 matched analogues with a probability > 0.7.
    :return: path of tab delimited text file with the selected motif, its fragments+probabilities and the MassQL query
    """
    file_path = Path(r"motif_massql_querries.txt")
    spectrum_file = open(file_path, "w")
    for index, row in df_motifs_to_frag.iterrows():
        if index in list_of_selected_motifs:
            spectrum_file.write("{0}    {1}    {2}".format(index, df_motifs_to_frag.at[index, "Fragment+Probability"],
                                                           make_MassQL_search(
                                                               df_motifs_to_frag.at[index, "Fragment+Probability"])))
            spectrum_file.write("\n")
    spectrum_file.close()
    return os.path.abspath(file_path)

def make_pdf(df_motifs_to_doc: pd.DataFrame, df_doc_to_smiles: pd.DataFrame, df_motifs_to_frag: pd.DataFrame,
                 threshold_amount_analogs=2) -> list:
    """selects certain motifs and makes an overview PDF with the Mass2Motif, fragments, and the associated analogs

    :param df_motifs_to_doc: Pandas dataframe, with the motifs as index and each associated document in a separate
    column. The last column is a list of lists of the associated document, the probability and the overlap score.
    :param df_doc_to_smiles: Pandas dataframe, with the associated document and their precursor m/z. Prediction
    certainty, precursor_mz_analog, and smiles of analogs are given for prediction certainties above threshold.
    :param df_motifs_to_frag: pandas dataframe, with the motifs as index and each feature and probability in a list in
    separate column. The last column contains a lists of list with each associated fragment and probability
    :return: list of motifs found in more than x (default=5) spectra in the dataframe and that contain at least one
    fragment/neutral loss with a prob > 0.05 and at least x (default=2) matched analogues with a probability > x
    (default=0.7).
    """

    pdf = FPDF()
    pdf.add_page()
    pdf.set_font("helvetica", size=10)
    list_of_selected_motifs = []
    for index, row in df_motifs_to_doc.iterrows():
        amount_of_smiles = 0
        # calculating the amount of smiles with a probability > 0.7 associated with each motif
        for cell in row:
            if cell in df_doc_to_smiles.index.values.tolist():
                if df_doc_to_smiles.at[cell, "smiles"] != "":
                    amount_of_smiles += 1
        # selection of Mass2Motifs that will be printed in PDF based on the amount of smiles with a probability > 0.7.
        # At least 2 or more Smiles are needed.
        if amount_of_smiles >= threshold_amount_analogs:
            pdf.multi_cell(200, 5, txt="\n", align='C')
            pdf.set_font("helvetica", "B", size=10)
            # writing the mass_2_motif name, which is the index of df_motifs_to_doc, in the file
            pdf.multi_cell(200, 5, txt="{0}\n".format(index), align='C')
            pdf.set_font("helvetica", size=10)
            # if the mass2motif is in the list of indexes of df_motifs_to_frag then look for the column with fragments
            # + probability of that mass2motif and sort that list of lists on the probability value
            if index in df_motifs_to_frag.index.values.tolist():
                list_of_selected_motifs.append(index)
                pdf.multi_cell(200, 5, txt="{0}\n".format(
                    sorted(df_motifs_to_frag.at[index, "Fragment+Probability"], key=lambda x: x[1], reverse=True)),
                               align='L')
                # https: // www.geeksforgeeks.org / python - sort - list - according - second - element - sublist /
                MassQL_query = make_MassQL_search(df_motifs_to_frag.at[index, "Fragment+Probability"])
                pdf.multi_cell(200, 5, txt="MassQL query: {0}\n".format(MassQL_query),
                               align='L')
            else:
                pdf.multi_cell(200, 5, txt="motif does not have fragments with probability > 0.05", align='L')
            # print the "Document+probability+Overlap" column for the selected index, which is the selected Mass2Motif
            pdf.set_font("helvetica", "I", size=10)
            pdf.multi_cell(200, 5,
                           txt="[[Spectrum number, probability, overlap score], [Spectrum number, probability, overlap score], etc.] :\n",
                           align='L')
            pdf.set_font("helvetica", size=10)
            pdf.multi_cell(200, 5, txt="{0}\n".format(df_motifs_to_doc.at[index, "Document+Probability+Overlap"]),
                           align='L')
            # print the amount of lists in the "Document+probability+Overlap" column for the selected index, which is
            # equal to the amount of documents the motif is associated with
            pdf.multi_cell(200, 5, txt="amount of spectra associated with motif: {0}\n".format(
                len(df_motifs_to_doc.at[index, "Document+Probability+Overlap"])), align='L')
            # each individual cell in df_motifs_to_doc is the document number associated with the index
            for cell in row:
                # if the cell value is an int, so if its really a spectrum number
                # (because there is also a list of lists in this dataframe)
                if isinstance(cell, int):
                    # if the spectrum number is in the index of df_doc_to_smiles
                    if cell in df_doc_to_smiles.index.values.tolist():
                        # if the smiles associated with the spectrum number is not empty, because most cells are empty
                        # because of a too low prediction probability of the analogue
                        if df_doc_to_smiles.at[cell, "smiles"] != "":
                            # then print the associated spectrum number that the smiles belongs to and the precursor m/z
                            # of the spectrum
                            pdf.multi_cell(200, 10, txt="spectrum number: {0}, precursor m/z: {1}\n".format(cell,
                                                                                                            df_doc_to_smiles.at[
                                                                                                                cell, "precursor_mz_query_spectrum"]),
                                           align='L')
                            # also use the function visualize_mol to print a picture of the structure of the analog
                            # using the smiles that is associated with the spectrum
                            pdf.image(visualize_mol(df_doc_to_smiles.at[cell, "smiles"]))
                            # print the prediction probability of the analog, and the mz of the smiles that is
                            # associated with the spectrum
                            pdf.multi_cell(200, 10, txt="analog m/z:{0}, prediction: {1}\n".format(
                                df_doc_to_smiles.at[cell, "precursor_mz_analog"],
                                df_doc_to_smiles.at[cell, "ms2query_model_prediction"]), align='L')
    pdf.output("output.pdf")
    return list_of_selected_motifs

def main() -> None:
    """Main function of this module"""

    #step 0: input of script
    path_to_file_with_MS2Query_csv= argv[1]
    path_file_with_MS2LDA_csv = argv[2]
    path_file_with_Motif_fragments_csv = argv[3]
    # step 1: rearrange MS2LDA file with mass2motifs and documents in dataframe
    df_motifs_to_doc=convert_Mass2Motifs_to_df(path_file_with_MS2LDA_csv, threshold_spectra_in_motif=5)
    #step 2: rearrange MS2Query file with documents and smiles of analogs in dataframe
    df_doc_to_smiles=convert_ms2query_to_df(path_to_file_with_MS2Query_csv, threshold_ms2query_pred=0.7)
    # step 3: rearrange MS2LDA fragment file with mass2motifs and fragments in dataframe
    df_motifs_to_frag=convert_fragments_in_motif_to_df(path_file_with_Motif_fragments_csv)
    #step 4: make a PDF where you can see each Mass2Motif, its fragments, and the associated analog structures
    list_of_selected_motifs=make_pdf(df_motifs_to_doc,df_doc_to_smiles,df_motifs_to_frag, threshold_amount_analogs=2)
    # step 5: make a table with the massql querries and their motifs and fragments
    print(make_file_with_massql_querries(df_motifs_to_frag, list_of_selected_motifs))

if __name__ == "__main__":
    main()