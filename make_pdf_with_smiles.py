#!/usr/bin/env python3
"""
Author: Anna-Marie Seelen
Studentnumber:1008970
Description: takes result csv files of MS2LDA and MS2Query and outputs a pdf that gives an overview of
motifs and analogs. In addition, it outputs a text file with the selected mass2motifs and their MassQL querries.
Usage: python3 *name_of_script* *path_to_file_with_MS2Query_csv* *path_file_with_MS2LDA_csv*
*path_file_with_Motif_fragments_csv* *path_with_mgf_file_from_GNPS_with_your_data*

    name_of_script: make_pdf_with_smiles.py
    path_to_file_with_MS2Query_csv: output file from the MS2Query.py script with document number and
    information about the matched compound
    path_file_with_MS2LDA_csv: file downloaded from MS2LDA.org which contains document number and name of the motif
    associated with the document
    path_file_with_Motif_fragments_csv: file downloaded from MS2LDA.org with contains the fragments and neutral losses
    and which motif these are associated with
    path_with_mgf_file_from_GNPS_with_your_data: file downloaded from the GNPS website after a molecular networking job.
    This file contains all the spectra from the inputted MzML format in a mgf style format.
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
from matchms.importing import load_from_mgf
from matchms.filtering import add_losses
from decimal import *

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
    df_motif_list_of_lists_feature = df_Motif_fragments[["Motif", "Fragment+Probability"]].groupby("Motif",
                                                                                                   as_index=True).aggregate(
        {"Fragment+Probability": list})
    # Make a second dataframe with the Motifs as index and a second column with a list containing the associated
    # features and probability in a separate column
    f = lambda x: 'fragment_{}+probability'.format(x + 1)
    df_motif_and_sep_fragments = df_motif_list_of_lists_feature["Fragment+Probability"].apply(pd.Series).fillna(
        '').rename(columns=f)
    # https://stackoverflow.com/questions/45107523/pandas-convert-list-of-lists-to-multiple-columns
    df_new=pd.concat([df_motif_and_sep_fragments, df_motif_list_of_lists_feature], axis=1)
    return df_new

def make_list_of_selected_motifs(df_motifs_to_doc: pd.DataFrame, df_doc_to_smiles: pd.DataFrame, df_motifs_to_frag: pd.DataFrame,
                 threshold_amount_analogs=2) -> list:
    """
    Generates a list containing motifs for which at least 2 with a probability > x (default=0.7) analogues were matched.

    :param df_motifs_to_doc: Pandas dataframe, with the motifs as index and each associated document in a separate
    column. The last column is a list of lists of the associated document, the probability and the overlap score.
    :param df_doc_to_smiles: Pandas dataframe, with the associated document and their precursor m/z. Prediction
    certainty, precursor_mz_analog, and smiles of analogs are given for prediction certainties above threshold.
    :param df_motifs_to_frag: pandas dataframe with the motifs as index and each associated document in a separate column. The last
    column is a list of lists of the associated document, the probability and the overlap score
    :param threshold_amount_analogs: int, the amount of matched analogues a motif needs to have for it to be deemed
    interesting (default=2).
    :return: list of motifs found in more than x (default=5) spectra in the dataframe and that contain at least one
    fragment/neutral loss with a prob > 0.05 and at least x (default=2) matched analogues with a probability > x
    (default=0.7).
    """
    list_of_selected_motifs = []
    for index, row in df_motifs_to_doc.iterrows():
        amount_of_smiles = 0
        # calculating the amount of smiles with a probability > 0.7 (default) associated with each motif
        for cell in row:
            if cell in df_doc_to_smiles.index.values.tolist():
                if df_doc_to_smiles.at[cell, "smiles"] != "":
                    amount_of_smiles += 1
        # selection of Mass2Motifs based on the amount of smiles with a probability > 0.7 (default).
        # At least 2 or more Smiles are needed for a mass2motif to be added to the list so that the validation can happen
        if amount_of_smiles >= threshold_amount_analogs:
            # if the mass2motif is in the list of indexes of df_motifs_to_frag then look for the column with fragments
            if index in df_motifs_to_frag.index.values.tolist():
                list_of_selected_motifs.append(index)
    return list_of_selected_motifs

def calculate_doc_ratio_for_feature(mgf_file, list_of_selected_motifs: list, df_motifs_to_frag,
                                     df_motifs_to_doc: pd.DataFrame, minimum_ratio=0) -> pd.DataFrame:
    """
    Calculates the ratio of the associated documents with each feature of each motif and selects features based on the ratio.

    :param mgf_file: str, path to file downloaded from the GNPS website after a molecular networking job.
    This file contains all the spectra from the inputted MzML format in a mgf style format.
    :param list_of_selected_motifs: list, list of motifs found in more than x (default=5) spectra in the dataframe and
    that contain at least one  fragment/neutral loss with a prob > 0.05 and at least x (default=2) matched analogues
    with a probability > x (default=0.7).
    :param df_motifs_to_frag: pandas dataframe with the motifs as index and each feature and probability in a list in
    separate column.The last column contains a lists of list with each associated fragment and probability
    :param df_motifs_to_doc: Pandas dataframe, with the motifs as index and each associated document in a separate
    column. The last column is a list of lists of the associated document, the probability and the overlap score.
    :param minimum_ratio: the minimum ratio of associated documents with feature/total associated documents with the
    motif that the feature with the highest ratio has to have for the mass2motif to be deemed interesting.
    :return: Pandas Dataframe with the selected motifs as an index and in one column a list of lists containing the
    feature, probability of the feature, the ratio of associated document with the feature, a list of the documents that
    contain the feature for every selected feature.
    """
    spectra = list(load_from_mgf(mgf_file))
    #create list with 0 to fill the pandas column for the features for now
    #convert 0 to string, because you get a Value error at the end of the function from pandas otherwise
    empty_list = ["0"]*(len(list_of_selected_motifs))
    # create pandas dataframe with the motifs that have 2 (default) smiles with prediction probability above 0.7 (default)
    data={"motif": list_of_selected_motifs, "Fragment+Probability+Ratio+Doc": empty_list}
    df_selected_motif_and_ratio = pd.DataFrame(data)
    df_selected_motif_and_ratio.set_index("motif", inplace=True)
    for motif in list_of_selected_motifs:
        features_list_of_lists_with_counts=[]
        for feature in df_motifs_to_frag.at[motif, "Fragment+Probability"]:
            documents_that_contain_feature=[]
            num_of_ass_doc_with_feature=0
            total_documents_num = len(df_motifs_to_doc.at[motif, "Document+Probability+Overlap"])
            #for every document associated with a motif we will see if the spectrum of the document contains the feature
            for document in df_motifs_to_doc.at[motif, "Document+Probability+Overlap"]:
                # for every document associated with each feature for each motif you want to check if the document
                # contains the feature
                for spectrum in spectra:
                    # look through all the matchms spectra select the spectrum that has the current document number
                    if int(spectrum.get("scans")) == int(document[0]):
                        # if the feature we are looking at is a loss
                        if re.search(r'loss', feature[0]) != None:
                            # first calculate all the losses in de spectrum based on the parent mass
                            spectrum=add_losses(spectrum)
                            for loss in range(len(spectrum.losses.mz)):
                                # get the losses with 2 decimal points and round down the number
                                rounded_loss = Decimal(spectrum.losses.mz[loss]).quantize(Decimal('.01'),
                                                                                             rounding=ROUND_DOWN)
                                # if a loss in the spectrum contains the loss in the feature then
                                # the document contains the loss!
                                if float(rounded_loss)==float(re.search(r'\_(.*\..{2}).*', feature[0]).group(1)):
                                    num_of_ass_doc_with_feature+=1
                                    documents_that_contain_feature.append(int(document[0]))
                        #if the feature we are looking at is a fragment
                        else:
                            for fragment in range(len(spectrum.peaks.mz)):
                                # get the fragment of the spectrum with 2 decimal points and round down the number
                                rounded_fragment = Decimal(spectrum.peaks.mz[fragment]).quantize(Decimal('.01'),
                                                                                                 rounding=ROUND_DOWN)
                                # if a fragment in the spectrum is the same as a fragment in the feature
                                # then the document contains the fragment!
                                if float(rounded_fragment) == float(re.search(r'\_(.*\..{2}).*', feature[0]).group(1)):
                                    num_of_ass_doc_with_feature+=1
                                    documents_that_contain_feature.append(int(document[0]))
            # calculate the ratio of the associated document with the feature
            ratio_of_ass_doc_with_feature=round((num_of_ass_doc_with_feature/total_documents_num),2)
            feature.append(ratio_of_ass_doc_with_feature)
            feature.append(documents_that_contain_feature)
            features_list_of_lists_with_counts.append(feature)

        #select the features that will be in the massql search based on the ratio of the associated doc with the feature
        #sort the list of lists of features from high to low ratio
        features_list_of_lists_with_counts = sorted(features_list_of_lists_with_counts, key=lambda x: x[2],
                                                    reverse=True)
        #if the highest feature ratio of a motif is below the minimum_ratio (default: 0)
        # the whole motif will be discarded from the dataframe
        if features_list_of_lists_with_counts[0][2] < float(minimum_ratio):
            df_selected_motif_and_ratio=df_selected_motif_and_ratio.drop(motif)
        #if the highest feature ratio of a motif is above the minumum_ratio the motif will be in the dataframe with one
        # or two features:
        else:
            if len(features_list_of_lists_with_counts)>1:
                # the second highest feature will only be included in the massql search if the ratio of above 0.5
                if features_list_of_lists_with_counts[1][2]>=0.5:
                    df_selected_motif_and_ratio.at[
                        motif, "Fragment+Probability+Ratio+Doc"] = features_list_of_lists_with_counts[:2]
                else:
                    df_selected_motif_and_ratio.at[
                        motif, "Fragment+Probability+Ratio+Doc"] = features_list_of_lists_with_counts[:1]
            else:
                df_selected_motif_and_ratio.at[
                    motif, "Fragment+Probability+Ratio+Doc"] = features_list_of_lists_with_counts[:1]

    return df_selected_motif_and_ratio

def make_MassQL_search(fragments: list) -> str:
    """
    Takes a list of lists containing fragments and probabilities and returns corresponding MassQL query

    :param fragments: list of lists containing fragments and probabilities belonging to a Mass2Motif
    :return: str, the corresponding MassQL query
    """
    query="QUERY scaninfo(MS2DATA) WHERE POLARITY = Positive " #MS1DATA doesn't give any results, so MS2DATA it is
    for fragment in fragments:
        for string in fragment:
            if type(string)==str:
                if re.search(r'loss', string) != None:
                    query+="AND MS2NL = {0}:TOLERANCEMZ={1} ".format(re.search(r'\_(.*)', string).group(1), 0.01)
                else:
                    query+="AND MS2PROD = {0}:TOLERANCEMZ={1} ".format(re.search(r'\_(.*)', string).group(1), 0.01)
            else:
                pass
    return query

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

def make_pdf(df_motifs_to_doc: pd.DataFrame, df_doc_to_smiles: pd.DataFrame,
                 df_selected_motif_and_ratio: pd.DataFrame) -> None:
    """Makes an overview PDF with the Mass2Motif, fragments, and the associated analogs

    :param df_motifs_to_doc: Pandas dataframe, with the motifs as index and each associated document in a separate
    column. The last column is a list of lists of the associated document, the probability and the overlap score.
    :param df_doc_to_smiles: Pandas dataframe, with the associated document and their precursor m/z. Prediction
    certainty, precursor_mz_analog, and smiles of analogs are given for prediction certainties above threshold.
    :param df_selected_motif_and_ratio: Pandas Dataframe with the selected motifs as an index and in one column a list
    of lists containing the feature, probability of the feature, the ratio of associated document with the feature,
    a list of the documents that contain the feature for every selected feature.
    """
    pdf = FPDF()
    pdf.add_page()
    pdf.set_font("helvetica", size=10)
    list_of_selected_motifs=list(df_selected_motif_and_ratio.index.values)
    for index, row in df_motifs_to_doc.iterrows():
        # selection of Mass2Motifs that will be printed in PDF based on the amount of smiles with a probability > 0.7.
        # At least 2 or more Smiles are needed.
        if index in list_of_selected_motifs:
            pdf.multi_cell(200, 5, txt="\n", align='C')
            pdf.set_font("helvetica", "B", size=10)
            # writing the mass_2_motif name, which is the index of df_motifs_to_doc, in the file
            pdf.multi_cell(200, 5, txt="{0}\n".format(index), align='C')
            pdf.set_font("helvetica", size=10)
            # if the mass2motif is in the list of indexes of df_motifs_to_frag then look for the column with fragments
            # + probability and count of that mass2motif
            if index in df_selected_motif_and_ratio.index.values.tolist():
                pdf.set_font("helvetica", "I", size=10)
                pdf.multi_cell(200, 5,
                               txt="[[feature, probability, ratio_of_spectra_with_feature, spectra number with feature], [feature, probability, ratio_of_spectra_with_feature, spectra number with feature], etc.] :\n",
                               align='L')
                pdf.multi_cell(200, 5, txt="{0}\n".format(df_selected_motif_and_ratio.at[index, "Fragment+Probability+Ratio+Doc"]),
                               align='L')
                # https: // www.geeksforgeeks.org / python - sort - list - according - second - element - sublist /
                MassQL_query = make_MassQL_search(df_selected_motif_and_ratio.at[index, "Fragment+Probability+Ratio+Doc"])
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

def make_file_with_massql_querries(df_selected_motif_and_ratio: pd.DataFrame) -> str:
    """Outputs a tab delimited file with the selected motifs, their fragments and their MassQL querries

    :param df_selected_motif_and_ratio: Pandas Dataframe with the selected motifs as an index and in one column a list
    of lists containing the feature, probability of the feature, the ratio of associated document with the feature,
    a list of the documents that contain the feature for every selected feature.
    :return: path of tab delimited text file with the selected motif, its fragments+probabilities and the MassQL query
    """
    file_path = Path(r"motif_massql_querries.txt")
    spectrum_file = open(file_path, "w")
    for index, row in df_selected_motif_and_ratio.iterrows():
        spectrum_file.write("{0}    {1}    {2}".format(index, df_selected_motif_and_ratio.at[index, "Fragment+Probability+Ratio+Doc"],
                                                           make_MassQL_search(
                                                               df_selected_motif_and_ratio.at[index, "Fragment+Probability+Ratio+Doc"])))
        spectrum_file.write("\n")
    spectrum_file.close()
    return os.path.abspath(file_path)

def main() -> None:
    """Main function of this module"""

    #step 0: input of script
    path_to_file_with_MS2Query_csv= argv[1]
    path_file_with_MS2LDA_csv = argv[2]
    path_file_with_Motif_fragments_csv = argv[3]
    path_to_gnps_output_mgf_file = argv[4]
    # step 1: rearrange MS2LDA file with mass2motifs and documents in dataframe
    df_motifs_to_doc=convert_Mass2Motifs_to_df(path_file_with_MS2LDA_csv, threshold_spectra_in_motif=5)
    #step 2: rearrange MS2Query file with documents and smiles of analogs in dataframe
    df_doc_to_smiles=convert_ms2query_to_df(path_to_file_with_MS2Query_csv, threshold_ms2query_pred=0.7)
    # step 3: rearrange MS2LDA fragment file with mass2motifs and fragments in dataframe
    df_motifs_to_frag=convert_fragments_in_motif_to_df(path_file_with_Motif_fragments_csv)
    # step 4: make a list of the selected motifs with at least x (default 2) reliable MS2Query results per motif
    list_of_selected_motifs = make_list_of_selected_motifs(df_motifs_to_doc, df_doc_to_smiles, df_motifs_to_frag,
                                                           threshold_amount_analogs=2)
    # step 5: calculate the amount of associated documents for each feature relative to the total amount of
    # documents associated with the respective motif (so calculate the ratio)
    df_selected_motif_and_ratio=calculate_doc_ratio_for_feature(path_to_gnps_output_mgf_file, list_of_selected_motifs, df_motifs_to_frag,
                                 df_motifs_to_doc)
    #step 4: make a PDF where you can see each Mass2Motif, its fragments, and the associated analog structures
    make_pdf(df_motifs_to_doc,df_doc_to_smiles,df_selected_motif_and_ratio)
    # step 5: make a table with the selected motif and their corresponding massql querries
    print(make_file_with_massql_querries(df_selected_motif_and_ratio))

if __name__ == "__main__":
    main()