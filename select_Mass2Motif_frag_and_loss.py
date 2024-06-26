#!/usr/bin/env python3
"""
Author: Anna-Marie Seelen
Studentnumber:1008970
Description: takes result csv files of MS2LDA and MS2Query and outputs a pdf that gives an overview of the selected
motifs and analogs. In addition, it outputs a text file with the selected mass2motifs and their MassQL queries.
Usage: python3 *name_of_script* *path_to_file_with_MS2Query_results_csv* *path_file_with_MS2LDA_Motif_spectrum_ids_csv*
*path_file_with_Motif_fragments_csv* *path_with_mgf_file_from_GNPS_molecular_network* *dir_to_save_pic_of_structures*

    name_of_script: make_pdf_with_smiles.py
    path_to_file_with_MS2Query_results_csv: output file from the MS2Query.py script with spectrum_id number and
    information about the matched compound
    path_file_with_MS2LDA_Motif_spectrum_ids_csv: file downloaded from MS2LDA.org which contains spectrum_id number and
    name of the motif associated with the spectrum_id
    path_file_with_Motif_fragments_csv: file downloaded from MS2LDA.org which contains the fragments and neutral losses
    and which motif these are associated with
    path_with_mgf_file_from_GNPS_molecular_network: file downloaded from the GNPS website after a molecular networking
    job. This file contains all the spectra from the inputted MzML format in a mgf style format.
    dir_to_save_pic_of_structures: directory to save png files of ms2query analogues for experimental spectra
"""

# import statements
from sys import argv
import pandas as pd
from rdkit.Chem.Draw import MolToImage
from fpdf import FPDF
import re
from pathlib import Path
import os
from matchms.importing import load_from_mgf
from matchms.filtering import add_losses
from decimal import *
import time
from sys import argv
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import DrawingOptions, MolDrawOptions

# functions
def convert_Mass2Motifs_to_df(path_file_with_Motif_spectrum_ids_csv: str, threshold_spectra_in_motif=4) -> pd.DataFrame:
    """
    Takes csv file from MS2LDA.org and rearranges the motifs and spectrum_ids and returns them in a pandas dataframe.

    :param path_file_with_Motif_spectrum_ids_csv: str, path to location of file with csv from MS2LDA containing motifs
    and identifiers of spectra associated with the motif
    :param threshold_spectra_in_motif: int, default=4, minimum amount of spectra associated with a motif for the motif
    to be included in the output dataframe
    :return: pandas dataframe with the motifs as index and each associated spectrum_id in a separate column. The last
    column is a list of lists with the associated spectrum_id, the probability and the overlap score
    """

    df_MS2LDA = pd.read_csv(path_file_with_Motif_spectrum_ids_csv, header=0)
    #add a new column the dataframe with a list of the spectrum_id, probability and overlap score for each spectrum_id
    df_MS2LDA.rename(columns = {'Document':'spectrum_id'}, inplace = True)
    df_MS2LDA["spectrum_id+Probability+Overlap"] = df_MS2LDA[
        ["spectrum_id", "Probability", "Overlap Score"]].values.tolist()
    # First make a dataframe with the Motifs as index and a second column with a list of lists containing the associated
    # spectrum_ids, probability and overlap scores
    df_motif_and_list = df_MS2LDA[["Motif", "spectrum_id+Probability+Overlap"]].groupby("Motif", as_index=True).aggregate(
        {"spectrum_id+Probability+Overlap": list})
    # We only want the motifs found in more than x spectra in the dataframe
    df_motif_and_list = df_motif_and_list[
        df_motif_and_list["spectrum_id+Probability+Overlap"].str.len() >= threshold_spectra_in_motif]
    # Make a second dataframe with the Motifs as index and a second column with a list containing the associated
    # spectrum_ids
    df_motif_and_spectralist = df_MS2LDA[["Motif", "spectrum_id"]].groupby("Motif", as_index=True).aggregate(
        {"spectrum_id": list})
    # Again we only want the motifs found in more than x spectra in the dataframe
    df_motif_and_spectralist = df_motif_and_spectralist[
        df_motif_and_spectralist['spectrum_id'].str.len() >= threshold_spectra_in_motif]
    # solution from: https://stackoverflow.com/questions/58297277/how-to-get-a-length-of-lists-in-pandas-dataframe
    # now we want each spectrum_id to be in a separate column in the dataframe
    f = lambda x: 'spectrum_id_{}'.format(x + 1)
    df_motif_and_spectra = pd.DataFrame(df_motif_and_spectralist.spectrum_id.values.tolist(), df_motif_and_spectralist.index,
                                    dtype=object).fillna('').rename(columns=f)
    # solution from:
    # https://stackoverflow.com/questions/44663903/pandas-split-column-of-lists-of-unequal-length-into-multiple-columns
    df_motifs_and_spectra_and_list = pd.concat([df_motif_and_spectra, df_motif_and_list], axis=1)
    return df_motifs_and_spectra_and_list

def convert_ms2query_to_df(path_to_file_with_MS2Query_csv: str, threshold_ms2query_pred=0.7) -> pd.DataFrame:
    """
    Takes resulting csv file from MS2Query, selects the spectrum_ids and analogs, returns them in a pandas dataframe.

    :param path_to_file_with_MS2Query_csv: str, path to location of file with result csv from MS2Query containing
    spectrum_ids and their analogs
    :param threshold_ms2query_pred: int, default=0.7, minimal ms2query prediction for the analogue to be included in the
    output dataframe
    :return: pandas dataframe with the associated spectrum_id and their precursor m/z. Prediction certainty,
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
    Takes a csv file from MS2LDA.org and rearranges the motifs and fragments and returns them in a pandas dataframe.

    :param path_file_with_MS2LDA_csv_fragments: str, path to location of csv file from MS2LDA containing motifs and
    Mass2Motif fragments and losses
    :return: pandas dataframe with the motifs as index and each fragment or loss and probability in a list in separate
    column. The last column contains a lists of list with each associated fragment and probability
    """

    df_Motif_fragments= pd.read_csv(path_file_with_MS2LDA_csv_fragments, header=0)
    # First make a column for each motif with a list of the probability and fragment or loss combined
    df_Motif_fragments.rename(columns={'Feature': 'frag_or_loss'}, inplace=True)
    df_Motif_fragments["Fragment+Probability"]=df_Motif_fragments[["frag_or_loss", "Probability"]].values.tolist()
    # Then make Motifs index and the second column with a list of lists containing the fragments and probability
    df_motif_list_of_lists_frag_or_loss = df_Motif_fragments[["Motif", "Fragment+Probability"]].groupby("Motif",
                                                                                                   as_index=True).aggregate(
        {"Fragment+Probability": list})
    # Make a second dataframe with the Motifs as index and a second column with a list containing the associated
    # mass2motif fragments and probability in a separate column
    f = lambda x: 'fragment_{}+probability'.format(x + 1)
    df_motif_and_sep_fragments = df_motif_list_of_lists_frag_or_loss["Fragment+Probability"].apply(pd.Series).fillna(
        '').rename(columns=f)
    # https://stackoverflow.com/questions/45107523/pandas-convert-list-of-lists-to-multiple-columns
    # concatenate the two dataframes
    df_motif_frag_separate_and_list=pd.concat([df_motif_and_sep_fragments, df_motif_list_of_lists_frag_or_loss], axis=1)
    return df_motif_frag_separate_and_list

def make_list_of_selected_motifs(df_motifs_to_spectra: pd.DataFrame, df_spectra_to_smiles: pd.DataFrame,
                                     df_motifs_to_frag: pd.DataFrame,
                                     threshold_amount_analogs=1) -> list:
    """
    Generates a list containing motifs for which at least 1 with a probability > x (default=0.7) analogues were matched.

    :param df_motifs_to_spectra: Pandas dataframe, with the motifs as index and each associated spectrum_id in a separate
    column. The last column is a list of lists with the associated spectrum_id, the probability and the overlap score
    :param df_spectra_to_smiles: Pandas dataframe, with the associated spectrum_id and their precursor m/z. Prediction
    certainty, precursor_mz_analog, and smiles of analogs are given for prediction certainties above threshold.
    :param df_motifs_to_frag: pandas dataframe with the motifs as index and each associated spectrum_id in a separate
    column. The last column is a list of lists of the associated spectrum_id, the probability and the overlap score
    :param threshold_amount_analogs: int, the amount of matched analogues a motif needs to have for it to be deemed
    interesting (default=1).
    :return: list of motifs found in more than x (default=4) spectra in the dataframe and that contain at least one
    fragment/neutral loss with a prob > 0.05 and at least x (default=1) matched analogues with a probability > x
    (default=0.7).

    A analog match was only required for Mass2Motifs that were unannotated in MotifDB!
    """
    list_of_selected_motifs = []
    for index, row in df_motifs_to_spectra.iterrows():
        amount_of_smiles = 0
        # calculating the amount of smiles with a probability > 0.7 (default) associated with each motif
        # if the index starts with motif the motif is unannotated and requires 1 spectrum (default) analogue
        if index.startswith("motif"):
            for cell in row:
                if cell in df_spectra_to_smiles.index.values.tolist():
                    if df_spectra_to_smiles.at[cell, "smiles"] != "":
                        amount_of_smiles += 1
            # selection of Mass2Motifs based on the amount of smiles with a probability > 0.7 (default).
            # At least 1 or more Smiles are needed for a mass2motif to be added to the list so that the validation can
            # happen
            if amount_of_smiles >= threshold_amount_analogs:
                # if the mass2motif is in the list of indexes of df_motifs_to_frag then look for the column with
                # fragments
                if index in df_motifs_to_frag.index.values.tolist():
                    list_of_selected_motifs.append(index)
        # else if the motif is annotated in MotifDB then it is selected straight away (also without a MS2Query match)
        else:
            list_of_selected_motifs.append(index)
    print(f"based on the selection criteria, the following motifs are selected: {list_of_selected_motifs}")
    return list_of_selected_motifs

def calculate_spectra_ratio_for_frag_or_loss(M2M_fragment_or_loss: list, motif: str,
                                                 M2M_frag_list_of_lists_with_counts: list, mgf_file: str,
                                                 df_motifs_to_spectra: pd.DataFrame) -> list:
    """
    Retrieves the associated spectrum_ids with a frag_or_loss and calculates the spectrum_id ratio for the frag_or_loss.

    :param M2M_fragment_or_loss: list, the M2M fragments and/or losses for which the spectrum_id ratio should be
    calculated
    :param motif: str, the motif name of the fragments and/or losses for which the spectrum_id ratio should be
    calculated
    :param M2M_frag_list_of_lists_with_counts: list, a list of lists with each list containing a fragment or loss
    associated with the motif, probability of the fragment or loss, the ratio of associated spectrum_id with the
    fragment or loss, and a list of the spectrum_ids that contain the fragment or loss.
    :param mgf_file: str, path to file downloaded from the GNPS website after a molecular networking job.
    This file contains all the spectra from the inputted MzML format in a mgf style format.
    :param df_motifs_to_spectra: Pandas dataframe, with the motifs as index and each associated spectrum_id in a
    separate column. The last column is a list of lists of the associated spectrum_id, the probability and the overlap
    score.
    :return: a list of lists with each list containing a fragment or loss associated with the motif, probability of the
    fragment or loss, the ratio of associated spectrum_id with the fragment or loss, a list of the spectrum_ids that
    contain the fragment or loss.
    """

    number_spectrum_ids_ass_with_motif = len(df_motifs_to_spectra.at[motif, "spectrum_id+Probability+Overlap"])
    spectra = list(load_from_mgf(mgf_file))
    spectrum_ids_that_contain_fragment = []
    num_of_ass_spectra_with_frag_or_loss = 0
    # for every spectrum_id associated with a motif we will see if the spectrum of the spectrum_id contains the fragment
    # or loss
    for spectrum_id in df_motifs_to_spectra.at[motif, "spectrum_id+Probability+Overlap"]:
        # for every spectrum_id associated with each fragment or loss for each motif you want to check if the spectrum_id
        # contains the fragment or loss
        for spectrum in spectra:
            # look through all the matchms spectra select the spectrum that has the current spectrum_id number
            if int(spectrum.get("scans")) == int(spectrum_id[0]):
                # if the fragment or loss we are looking at is a loss
                if re.search(r'loss', M2M_fragment_or_loss[0]) != None:
                    # first calculate all the losses in de spectrum based on the parent mass
                    spectrum = add_losses(spectrum)
                    for loss in range(len(spectrum.losses.mz)):
                        # get the losses with 2 decimal points and round down the number
                        rounded_loss = Decimal(spectrum.losses.mz[loss]).quantize(Decimal('.01'),
                                                                                  rounding=ROUND_DOWN)
                        # if a loss in the spectrum contains the loss in the fragment or loss then
                        # the spectrum_id contains the loss!
                        if float(rounded_loss) == float(re.search(r'\_(.*\..{2}).*', M2M_fragment_or_loss[0]).group(1)):
                            num_of_ass_spectra_with_frag_or_loss += 1
                            spectrum_ids_that_contain_fragment.append(int(spectrum_id[0]))
                # if the fragment or loss we are looking at is a fragment
                else:
                    for fragment in range(len(spectrum.peaks.mz)):
                        # get the fragment of the spectrum with 2 decimal points and round down the number
                        rounded_fragment = Decimal(spectrum.peaks.mz[fragment]).quantize(Decimal('.01'),
                                                                                         rounding=ROUND_DOWN)
                        # if a fragment in the spectrum is the same as a fragment in the fragment or loss
                        # then the spectrum_id contains the fragment!
                        if float(rounded_fragment) == float(
                                re.search(r'\_(.*\..{2}).*', M2M_fragment_or_loss[0]).group(1)):
                            num_of_ass_spectra_with_frag_or_loss += 1
                            spectrum_ids_that_contain_fragment.append(int(spectrum_id[0]))
    # calculate the ratio of the associated spectrum_id with the Mass2Motif fragment or loss
    ratio_of_ass_spectra_with_frag_or_loss = round(
        (num_of_ass_spectra_with_frag_or_loss / number_spectrum_ids_ass_with_motif), 2)
    M2M_fragment_or_loss.append(ratio_of_ass_spectra_with_frag_or_loss)
    M2M_fragment_or_loss.append(spectrum_ids_that_contain_fragment)
    M2M_frag_list_of_lists_with_counts.append(M2M_fragment_or_loss)
    return M2M_frag_list_of_lists_with_counts

def select_motif_frag_based_on_spectra_ratio(mgf_file: str, list_of_selected_motifs: list,
                                                 df_motifs_to_frag: pd.DataFrame,
                                                 df_motifs_to_spectra: pd.DataFrame, minimum_ratio=0) -> pd.DataFrame:
    """
    Calculates the ratio of the associated spectrum_ids with each fragment or loss of each motif and selects fragments
    and/or losses based on the ratio.

    :param mgf_file: str, path to file downloaded from the GNPS website after a molecular networking job.
    This file contains all the spectra from the inputted MzML format in a mgf style format.
    :param list_of_selected_motifs: list, list of motifs found in more than x (default=4) spectra in the dataframe and
    that contain at least one  fragment/neutral loss with a prob > 0.05 and at least x (default=1) matched analogues
    with a probability > 0.05 (default=0.7).
    :param df_motifs_to_frag: pandas dataframe with the motifs as index and each fragment or loss and probability in a
    list in separate column.The last column contains a lists of list with each associated fragment and probability
    :param df_motifs_to_spectra: Pandas dataframe, with the motifs as index and each associated spectrum_id in a separate
    column. The last column is a list of lists of the associated spectrum_id, the probability and the overlap score.
    :param minimum_ratio: the minimum ratio of associated spectrum_ids with fragment or loss/total associated
    spectrum_ids with the motif that the fragment or loss with the highest ratio has to have for the mass2motif to be
    deemed interesting.
    :return: Pandas Dataframe with the selected motifs as an index and in one column a list of lists containing the
    fragment or loss, probability of the fragment or loss, the ratio of associated spectrum_id with the fragment or loss
    , a list of the spectrum_ids that contain the fragment or loss for every selected fragment or loss.
    """

    #create a pandas results column
    #create list with 0 to fill the pandas column for the fragments and/or losses for now
    #convert 0 to string, because you get a Value error at the end of the function from pandas otherwise
    empty_list = ["0"]*(len(list_of_selected_motifs))
    empty_list_1 = ["0"] * (len(list_of_selected_motifs))
    # create pandas dataframe with the selected motifs
    data = {"motif": list_of_selected_motifs, "Fragment+Probability+Ratio+spectra": empty_list,
            "All_Fragment+Probability+Ratio+spectra": empty_list_1}
    df_selected_motif_and_ratio = pd.DataFrame(data)
    df_selected_motif_and_ratio.set_index("motif", inplace=True)
    for motif in list_of_selected_motifs:
        # list with all fragments and/or losses that the motif contains
        frag_loss_list_of_lists_with_counts=[]
        # for every fragment or loss you want to calculate the ratio of associated spectrum_ids with the fragment or
        # loss with respect to the amount of spectrum_ids associated with the whole motif
        for frag_or_loss in df_motifs_to_frag.at[motif, "Fragment+Probability"]:
            calculate_spectra_ratio_for_frag_or_loss(frag_or_loss, motif, frag_loss_list_of_lists_with_counts, mgf_file,
                                            df_motifs_to_spectra)
        #select the fragments and/or losses that will be in the massql search based on the ratio of the associated
        # spectra with the frag_or_loss sort the list of lists of fragments and/or losses from high to low ratio
        frag_loss_list_of_lists_with_counts = sorted(frag_loss_list_of_lists_with_counts, key=lambda x: x[2],
                                                    reverse=True)
        df_selected_motif_and_ratio.at[
            motif, "All_Fragment+Probability+Ratio+spectra"] = frag_loss_list_of_lists_with_counts
        #if the highest fragment or loss ratio of a motif is below the minimum_ratio (default: 0)
        #the whole motif will be discarded from the dataframe
        if frag_loss_list_of_lists_with_counts[0][2] < float(minimum_ratio):
            df_selected_motif_and_ratio=df_selected_motif_and_ratio.drop(motif)
        #if the highest fragment or loss ratio of a motif is above the minumum_ratio the motif will be in the dataframe
        # with one or two fragments and/or losses:
        else:
            if len(frag_loss_list_of_lists_with_counts)>1:
                # the second highest fragment or loss will only be included in the massql search if the ratio of above 0.5
                if frag_loss_list_of_lists_with_counts[1][2]>=0.5:
                    df_selected_motif_and_ratio.at[
                        motif, "Fragment+Probability+Ratio+spectra"] = frag_loss_list_of_lists_with_counts[:2]
                else:
                    df_selected_motif_and_ratio.at[
                        motif, "Fragment+Probability+Ratio+spectra"] = frag_loss_list_of_lists_with_counts[:1]
            else:
                df_selected_motif_and_ratio.at[
                    motif, "Fragment+Probability+Ratio+spectra"] = frag_loss_list_of_lists_with_counts[:1]
    return df_selected_motif_and_ratio

def cal_amount_of_spectrum_ass_with_amount_of_frag(df_selected_motif_and_ratio: pd.DataFrame) -> None:
    """
    calculates how many spectra in the selected mass2motifs are associated with 1 fragment_or_loss and with more

    :param df_selected_motif_and_ratio: Pandas Dataframe with the selected motifs as an index and in one column a list
    of lists containing the fragment or loss, probability of the fragment or loss, the ratio of associated spectrum_id
    with the fragment or loss, a list of the spectrum_ids that contain the fragment or loss for every selected fragment
    or loss.
    :return: None

    This is done before the selection of Mass2Motif fragments based on counts, so all the Mass2Motif fragments of the
    selected Mass2Motifs are included in the analysis
    """

    list_all_spectrum_ids_in_all_M2M = []
    for index, row in df_selected_motif_and_ratio.iterrows():
        print(index)
        frag_loss_list_of_lists_with_counts = row['All_Fragment+Probability+Ratio+spectra']
        #calculate the amount of spectra that are associated with more than 1 fragment or loss within the same motif
        list_with_all_spectrum_in_M2M = []
        for frag_or_loss_list in frag_loss_list_of_lists_with_counts:
            associated_spectra_list = frag_or_loss_list[3]
            for num in associated_spectra_list:
                list_with_all_spectrum_in_M2M.append(num)
        for spectra in set(list_with_all_spectrum_in_M2M):
            list_all_spectrum_ids_in_all_M2M.append(spectra)
        set_multiple = set()
        count_single = 0
        for i in list_with_all_spectrum_in_M2M:
            # if the spectrum identifier is multiple times in the list of spectra associated with the mass2motif
            # the spectrum is associated with more than 1 mass2motif fragment or loss
            if [spectra for spectra in list_with_all_spectrum_in_M2M].count(i) > 1:
                set_multiple.add(i)
            # else the spectrum is only associated with one mass2motif fragment or loss
            else:
                count_single += 1
        count_multiple = len(set_multiple)
        list_with_counts = [count_single, count_multiple]
        print(f"count_single, count_multiple: {list_with_counts}")

    # determine which spectra are associated with more than 2 Mass2Motifs
    spectrum_ids_in_more_than_1_motif = []
    for i in list_all_spectrum_ids_in_all_M2M:
        if [i for i in list_all_spectrum_ids_in_all_M2M].count(i) > 1:
            spectrum_ids_in_more_than_1_motif.append(i)
    print(f"spectrum_ids that are ass with more than 1 selected motif {spectrum_ids_in_more_than_1_motif}")

    return None

def make_MassQL_search(fragments: list, intensity_threshold=0.8, mz_tolerance=0.01) -> str:
    """
    Takes a list of lists with selected Mass2Motif fragments and probabilities and returns corresponding MassQL query

    :param fragments: list of lists containing fragments and probabilities belonging to a Mass2Motif
    :param intensity_threshold: float, the relative intensity the selected fragments have to be in a spectrum to be
    matched (default: 0.8)
    :param mz_tolerance: float, the deviation in Da from the m/z of the Mass2Motif fragment or loss that MassQL will
    allow for a match in a spectrum (default: 0.01)
    :return: str, the corresponding MassQL query
    """
    # the beginning of the query
    query="QUERY scaninfo(MS2DATA) WHERE POLARITY = Positive " #MS1DATA doesn't give any results, so MS2DATA it is
    for fragment in fragments:
        for string in fragment:
            if type(string)==str:
                if re.search(r'loss', string) != None:
                    query += "AND MS2NL = {0}:TOLERANCEMZ={1}:INTENSITYVALUE={2} ".format(
                        re.search(r'\_(.*)', string).group(1), mz_tolerance, intensity_threshold)
                else:
                    query += "AND MS2PROD = {0}:TOLERANCEMZ={1}:INTENSITYVALUE={2} ".format(
                        re.search(r'\_(.*)', string).group(1), mz_tolerance, intensity_threshold)
            else:
                pass
    return query

def visualize_mol(smiles: str):
    """
    Takes a smiles as input and outputs an PIL<PNG> image

    :param smiles: str, a smiles string that corresponds to a molecular structure
    :return: PIL<PNG> image, image of structure corresponding with the smiles which can be inserted in a PDF

    function adapted from: https://pchanda.github.io/See-substructure-in-molecule/
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    Chem.Kekulize(mol)
    img = MolToImage(mol, size=(200, 200), fitImage=True)
    return img

def make_pdf(df_motifs_to_spectra: pd.DataFrame, df_spectra_to_smiles: pd.DataFrame,
                 df_selected_motif_and_ratio: pd.DataFrame, path_to_save_png_of_structures) -> None:
    """Makes an overview PDF with the Mass2Motif, fragments, and the associated analogs

    :param df_motifs_to_spectra: Pandas dataframe, with the motifs as index and each associated spectrum_id in a
    separate column. The last column is a list of lists of the associated spectrum_id, the probability and the overlap
    score.
    :param df_spectra_to_smiles: Pandas dataframe, with the associated spectrum_id and their precursor m/z. Prediction
    certainty, precursor_mz_analog, and smiles of analogs are given for prediction certainties above threshold.
    :param df_selected_motif_and_ratio: Pandas Dataframe with the selected motifs as an index and in one column a list
    of lists containing the fragment or loss, probability of the fragment or loss, the ratio of associated spectrum_id
    with the fragment or loss, a list of the spectrum_ids that contain the fragment or loss for every selected fragment
    or loss.
    """
    pdf = FPDF()
    pdf.add_page()
    pdf.set_font("helvetica", size=10)
    list_of_selected_motifs=list(df_selected_motif_and_ratio.index.values)
    for index, row in df_motifs_to_spectra.iterrows():
        # selection of Mass2Motifs that will be printed in PDF based on the amount of smiles with a probability > 0.7.
        # At least 2 or more Smiles are needed.
        if index in list_of_selected_motifs:
            pdf.multi_cell(200, 5, txt="\n", align='C')
            pdf.set_font("helvetica", "B", size=10)
            # writing the mass_2_motif name, which is the index of df_motifs_to_spectra, in the file
            pdf.multi_cell(200, 5, txt="{0}\n".format(index), align='C')
            pdf.set_font("helvetica", size=10)
            # if the mass2motif is in the list of indexes of df_motifs_to_frag then look for the column with fragments
            # + probability and count of that mass2motif
            if index in df_selected_motif_and_ratio.index.values.tolist():
                pdf.set_font("helvetica", "I", size=10)
                pdf.multi_cell(200, 5,
                               txt="[[frag_or_loss, probability, ratio_of_spectra_with_frag_or_loss, spectra number with frag_or_loss], [frag_or_loss, probability, ratio_of_spectra_with_frag_or_loss, spectra number with frag_or_loss], etc.] :\n",
                               align='L')
                pdf.multi_cell(200, 5, txt="{0}\n".format(
                    df_selected_motif_and_ratio.at[index, "All_Fragment+Probability+Ratio+spectra"]),
                               align='L')
                pdf.multi_cell(200, 5,
                               txt="selected fragments and/or losses:\n",
                               align='L')
                pdf.multi_cell(200, 5, txt="{0}\n".format(df_selected_motif_and_ratio.at[index, "Fragment+Probability+Ratio+spectra"]),
                               align='L')
                # https: // www.geeksforgeeks.org / python - sort - list - according - second - element - sublist /
                MassQL_query = make_MassQL_search(df_selected_motif_and_ratio.at[index, "Fragment+Probability+Ratio+spectra"])
                pdf.multi_cell(200, 5, txt="MassQL query: {0}\n".format(MassQL_query),
                               align='L')
            else:
                pdf.multi_cell(200, 5, txt="motif does not have fragments with probability > 0.05", align='L')
            # print the "spectrum_id+probability+Overlap" column for the selected index, which is the selected Mass2Motif
            pdf.set_font("helvetica", "I", size=10)
            pdf.multi_cell(200, 5,
                           txt="[[Spectrum number, probability, overlap score], [Spectrum number, probability, overlap score], etc.] :\n",
                           align='L')
            pdf.set_font("helvetica", size=10)
            pdf.multi_cell(200, 5, txt="{0}\n".format(df_motifs_to_spectra.at[index, "spectrum_id+Probability+Overlap"]),
                           align='L')
            # print the amount of lists in the "spectrum_id+probability+Overlap" column for the selected index, which is
            # equal to the amount of spectrum_ids the motif is associated with
            pdf.multi_cell(200, 5, txt="amount of spectra associated with motif: {0}\n".format(
                len(df_motifs_to_spectra.at[index, "spectrum_id+Probability+Overlap"])), align='L')
            # each individual cell in df_motifs_to_spectra is the spectrum_id number associated with the index
            for cell in row:
                # if the cell value is an int, so if its really a spectrum number
                # (because there is also a list of lists in this dataframe)
                if isinstance(cell, int):
                    # if the spectrum number is in the index of df_spectra_to_smiles
                    if cell in df_spectra_to_smiles.index.values.tolist():
                        # if the smiles associated with the spectrum number is not empty, because most cells are empty
                        # because of a too low prediction probability of the analogue
                        if df_spectra_to_smiles.at[cell, "smiles"] != "":
                            # then print the associated spectrum number that the smiles belongs to and the precursor m/z
                            # of the spectrum
                            pdf.multi_cell(200, 10, txt="spectrum number: {0}, precursor m/z: {1}\n".format(cell,
                                                                                                            df_spectra_to_smiles.at[
                                                                                                                cell, "precursor_mz_query_spectrum"]),
                                           align='L')
                            # also use the function visualize_mol to print a picture of the structure of the analog
                            # using the smiles that is associated with the spectrum
                            pdf.image(visualize_mol(df_spectra_to_smiles.at[cell, "smiles"]))
                            # cell is spectrum number and the index is the name of the motif
                            vis_MS2Query_mol(path_to_save_png_of_structures, df_spectra_to_smiles.at[cell, "smiles"], cell, index)
                            # print the prediction probability of the analog, and the mz of the smiles that is
                            # associated with the spectrum
                            pdf.multi_cell(200, 10, txt="analog m/z:{0}, prediction: {1}\n".format(
                                df_spectra_to_smiles.at[cell, "precursor_mz_analog"],
                                df_spectra_to_smiles.at[cell, "ms2query_model_prediction"]), align='L')
    pdf.output("output.pdf")

def make_file_with_massql_querries(df_selected_motif_and_ratio: pd.DataFrame) -> str:
    """Outputs a tab delimited file with the selected motifs, their fragments and their MassQL querries

    :param df_selected_motif_and_ratio: Pandas Dataframe with the selected motifs as an index and in one column a list
    of lists containing the fragment or loss, probability of the fragment or loss, the ratio of associated spectrum_id
    with the fragment or loss, a list of the spectrum_ids that contain the fragment or loss for every selected fragment
    or loss.
    :return: path of tab delimited text file with the selected motif, its fragments+probabilities and the MassQL query
    """
    file_path = Path(r"motif_massql_querries.txt")
    motif_query_file = open(file_path, "w")
    for index, row in df_selected_motif_and_ratio.iterrows():
        motif_query_file.write("{0}\t{1}\t{2}".format(index, df_selected_motif_and_ratio.at[index, "Fragment+Probability+Ratio+spectra"],
                                                           make_MassQL_search(
                                                               df_selected_motif_and_ratio.at[index, "Fragment+Probability+Ratio+spectra"])))
        motif_query_file.write("\n")
    motif_query_file.close()
    return os.path.abspath(file_path)

def vis_MS2Query_mol(path_to_save_png: str, smiles: str, spectrum_num: int, motif: str) -> None:
    """
    Takes the SMILES string and makes a png file of the structure

    :param smiles: string, SMILES of an analogue of a query spectrum
    :param spectrum_num: int, the identifier of the experimental spectrum for which the analogue was identified
    :param motif: str, name of the motif to which the spectrum was associated
    :return: None
    """
    opts = MolDrawOptions()
    opts.updateAtomPalette({k: (0, 0, 0) for k in DrawingOptions.elemDict.keys()})

    mol = Chem.MolFromSmiles(smiles)
    Draw.MolToFile(mol, f"{path_to_save_png}/{spectrum_num}_{motif}.png", options=opts)
    return None

def main() -> None:
    """Main function of this module"""
    #step 0: input of script
    before_script = time.perf_counter()
    path_to_file_with_MS2Query_csv= argv[1]
    path_file_with_Motif_spectrum_ids_csv = argv[2]
    path_file_with_Motif_fragments_csv = argv[3]
    path_to_gnps_output_mgf_file = argv[4]
    path_to_save_png_of_structures = argv[5]
    # step 1: rearrange MS2LDA file with mass2motifs and spectrum_ids in dataframe
    df_motifs_to_spectra=convert_Mass2Motifs_to_df(path_file_with_Motif_spectrum_ids_csv, threshold_spectra_in_motif=4)
    #step 2: rearrange MS2Query file with spectrum_ids and smiles of analogs in dataframe
    df_spectra_to_smiles=convert_ms2query_to_df(path_to_file_with_MS2Query_csv, threshold_ms2query_pred=0.7)
    # step 3: rearrange MS2LDA fragment file with mass2motifs and fragments in dataframe
    df_motifs_to_frag=convert_fragments_in_motif_to_df(path_file_with_Motif_fragments_csv)
    # step 4: make a list of the selected motifs with at least x (default 2) reliable MS2Query results per motif
    list_of_selected_motifs = make_list_of_selected_motifs(df_motifs_to_spectra, df_spectra_to_smiles, df_motifs_to_frag,
                                                           threshold_amount_analogs=1)
    # step 5: calculate the amount of associated spectrum_ids for each fragment or loss relative to the total amount of
    # spectrum_ids associated with the respective motif (so calculate the ratio)
    df_selected_motif_and_ratio=select_motif_frag_based_on_spectra_ratio(path_to_gnps_output_mgf_file,
                                                                         list_of_selected_motifs, df_motifs_to_frag,
                                                                         df_motifs_to_spectra)
    # step 6: calculate the amount of spectra in the selected Mass2Motifs that are associated with only one mass2motif
    # fragment or loss
    cal_amount_of_spectrum_ass_with_amount_of_frag(df_selected_motif_and_ratio)
    #step 7: make a PDF where you can see each Mass2Motif, its fragments, and the associated analog structures
    make_pdf(df_motifs_to_spectra,df_spectra_to_smiles,df_selected_motif_and_ratio, path_to_save_png_of_structures)
    # step 8: make a table with the selected motif and their corresponding massql queries
    print(f"path to file with the massql queries: {make_file_with_massql_querries(df_selected_motif_and_ratio)}")
    after_script = time.perf_counter()
    print("how long the total script took {0}".format(after_script - before_script))

if __name__ == "__main__":
    main()
