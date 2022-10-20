import sys
#import ALL_GNPS_210409_positive_processed_annotated_CF_NPC_classes.txt
import ntpath
import pyarrow.feather as feather
import os
import pyteomics
import re
# import rdkit.Chem as Chem
# from rdkit.Chem.Draw import MolToImage
#from fpdf import FPDF
import json
from pandas import json_normalize
import pandas as pd
import ast
import subprocess
from pathlib import Path
import re
import os, glob
import matchms
from matchms import Scores, Spectrum
import json
from typing import List
import numpy
from matchms import Spectrum
import csv
from matchms.importing import load_from_mgf
from matchms.exporting import save_as_mgf
from typing import List


def annotate_peaks(spectrum_file, smiles):
    """
    Annotates the MS2 peaks given a smiles and a fragmentation spectrum
    :return:
    """
    #print(df.loc["CCMSLIB00006126912"])
    return None

def annotate_peaks(spectrum_file_name, smiles, identifier, abs_mass_tol=0.01):
    """
    Annotates the MS2 peaks given a smiles and a fragmentation spectrum
    :return:
    """
    # print(df.loc["CCMSLIB00004678842"])
    # cfm-annotate.exe <smiles_or_inchi> <spectrum_file> **<id>** <ppm_mass_tol> <abs_mass_tol> 0.01<param_file> <config_file> <output_file
    file_path_out = Path(r"/lustre/BIF/nobackup/seele006/cfm_annotation_out_{0}".format(identifier))
    out_fn = "cfm_annotation_out_{0}".format(identifier)
    if os.path.exists(out_fn):
        return out_fn
    cmd = 'cfm-annotate \'{0}\' {1} {2} -abs_mass_tol {3} -output_file {4}' \
            .format(smiles, spectrum_file_name, identifier, abs_mass_tol, out_fn)
    e = subprocess.check_call(cmd, shell=True)
    print("EXIT STATUS AND TYPE", e, type(e))
    print("hi")
    return out_fn

def visualize_mol(smiles: str) -> None:
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

def look_for(cfm_file,motif,fragments):
    """
    # First parse the cfm_file because you only want the middle part
    #Then for each fragment of the motif you want you want to select the smiles
    Selection of the annotated peaks that are relevant for the mass2motif
    #TODO: think about how you will do this for neutral loss...
    :return:
    """
    return None

def select_the_most_likely_option_for_one_motif():
    #look through all the possible smiles for a fragment and select the one that is there most often
    # return the motif with the with fragment, probability and most likely smiles
    return None

def make_spectrum_file_for_id_matchms(gnps_pickled_lib, spectrum_id, path_to_store_spectrum_files):
    #object=list(matchms.importing.load_scores.scores_from_pickle(gnps_pickled_lib))
    obj = pd.read_pickle(gnps_pickled_lib)
    # with open(gnps_pickled_lib, 'r') as f:
    #     spectra=pickletools.(f)
    # spectra=gnps_pickled_libS
    # print(spectra.peaks.mz[0])
    for spectrum in obj:
        if spectrum.get("spectrum_id") == spectrum_id:
            file_path = Path(r"{0}/spectrum_file_{1}_new.txt".format(path_to_store_spectrum_files, spectrum_id))
            print(os.path.abspath(file_path))
            spectrum_file=open(file_path, "w")
            for i in range(1):  #change to 3 for cfm-annotate
                #spectrum_file.write("energy{0}\n".format(i))
                for fragment in range(len(spectrum.peaks.mz)):
                    spectrum_file.write("{0} {1}".format(spectrum.peaks.mz[fragment], spectrum.peaks.intensities[fragment]))
                    spectrum_file.write("\n")
            spectrum_file.close()
    #object=list(matchms.importing.load_from_pickle)
    #print(object)
    return None

def make_spectrum_file_for_id(df_json, spectrum_id, path_to_store_spectrum_files):
    """Takes a nested sorted list and outputs a tab delimited file

    alignment_list: nested list with families and alignment lenghts
    return: tab delimited text file with the contents of each sub list on a line
    """
    list_of_lists=ast.literal_eval(df_json.loc[spectrum_id, "peaks_json"])
    file_path = Path(r"{0}/spectrum_file_{1}.txt".format(path_to_store_spectrum_files,spectrum_id))
    # if os.path.exists(file_path):
    #     print(os.path.abspath(file_path))
    #     return os.path.abspath(file_path)
    spectrum_file=open(file_path, "w")
    for i in range(3):
        spectrum_file.write("energy{0}\n".format(i))
        for sub_list in list_of_lists:
            #if sub_list[1]>600: #dit is ff een tussen oplossing om een file te krijgen waar je iets mee kan!
            spectrum_file.write("{0} {1}".format(sub_list[0], sub_list[1]))
            spectrum_file.write("\n")
    spectrum_file.close()
    return os.path.abspath(file_path)

def parse_input(filename):
    """Parses argonaut formatted file to extract accession number, organism name and DNA_seq and stores those in dict.

    filename: str, name of argonaut formatted input file
    return: nested dictionary with {accession_number:{organism_name:DNA-seq}}
    """

    lines=(open(filename))
    record_bool = False
    records=[]
    record = ""
    for line in lines:
        line=line.strip()
        line = line.replace('\n', '')
        line = line.replace('\t', '')
        if line.startswith("BEGIN IONS"):
            #line = line.replace("ACCESSION   ", "")
            record_bool=True
        elif line.startswith("END IONS"):
            record_bool=False
            records.append(record)
            record=""
        # elif line.startswith("ORGANISM"):
        #     organism_dict = {}
        #     line = line.replace("ORGANISM  ", "")
        #     key = list(gb_dict)[-1]
        #     organism_dict[line] = ""
        #     gb_dict[key] = organism_dict
        # elif "ORIGIN" in line:
        #     origin=True
        # elif "//" in line:
        #     origin=False
        if record_bool:
            record+=line
    #print(records)

    # dict={}
    for record in records:
        key = re.search(r'190.119(.*)520', record)
        if key is not None:
            print(record)
        # key=re.search(r'NAME=(.*)',record)
        # smiles=re.search(r'SMILES=(.*)',record)
        # mass=re.search(r'PEPMASS=(.*)',record)
        # if key is not None:
        #     dict[key.group(1)]=[]
        # if smiles is not None:
        #     last_key = list(dict)[-1]
        #     dict[last_key].append(smiles.group(1))
        # if mass is not None:
        #     last_key = list(dict)[-1]
        #     dict[last_key].append(mass.group(1))
        #INCHI=
    # print(dict)
            # line=line.replace(" ", "")
            # line=''.join(filter(lambda ch: not ch.isdigit(), line))
            # #https://www.studytonight.com/python-howtos/remove-numbers-from-string-in-python
            # line = line.replace("ORIGIN", "")
            # key = list(gb_dict)[-1]
            # dict=gb_dict[key]
            # last_key=list(dict)[-1]
            # gb_dict[key][last_key]+=line
    return None

def mgf_to_spectra(mgf_file):
    with open(mgf_file, 'r') as spectra_file:
        spectra_from_file = list(load_from_mgf(spectra_file))
        for spectrum in spectra_from_file:
            if spectrum.get("spectrum_id") == spectrum_id:

                save_as_mgf(spectrum, filename)

def save_json_as_csv(json_file):
    with open(json_file) as json_file:
        jsondata = json.load(json_file)

    data_file = open('/lustre/BIF/nobackup/seele006/MassQL_speclib/HMDB.csv', 'w', newline='')
    csv_writer = csv.writer(data_file)

    count = 0
    for data in jsondata:
        if count == 0:
            header = data.keys()
            print(header)
            csv_writer.writerow(header)
            count += 1
        csv_writer.writerow(data.values())

    data_file.close()
    return None

def save_as_json(spectrums: List[Spectrum], filename: str):
    """Save spectrum(s) as json file.
    :py:attr:`~matchms.Spectrum.losses` of spectrum will not be saved.
    Arguments:
    ----------
    spectrums:
        Expected input are match.Spectrum.Spectrum() objects.
    filename:
        Provide filename to save spectrum(s).
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

            json.dump(spectrum_dict, fout)
            if i < len(spectrums) - 1:
                fout.write(",")
        fout.write("]")

def make_json_file(pickle_file, path_to_store_json_file):
    # if os.path.exists(path_to_store_json_file):
    #      print(os.path.abspath(path_to_store_json_file))
    obj = pd.read_pickle(pickle_file)
    spectrum_list=[]
    for spectrum in obj:
        spectrum_list.append(spectrum)
    save_as_json(spectrum_list,path_to_store_json_file)
    return None

    with open(mgf_file, 'r') as spectra_file:
        spectra_from_file = list(load_from_mgf(spectra_file))
        for spectrum in spectra_from_file:
            #print(spectrum.get("spectrumid"))
            if spectrum.get("spectrumid") == spectrum_id:
                # in this function I get an error but don't know why
                save_as_mgf(spectrum, file_path)
    return None

def new_dataframe(df_massql_matches,df_json):
    """
    Makes a new dataframe based on the df_json and df_massql_matches with the smiles and the peaks for each match.

    :param df_massql_matches: returns a Pandas dataframe with the spectrum ids of the spectra as the index which contain
    the characteristics of the query.
    :param df_json: pandas Dataframe with one row for each spectrum and the information belonging to the spectrum in
    separate columns.
    :return:
    """

    df=pd.merge(df_massql_matches["precmz"],df_json[["Precursor_MZ","Smiles", "peaks_json"]],left_index=True, right_index=True)
    # We don't want a dataframe with the smiles anymore, because we are not using cfm-annotate. We just need the identifiers.
    # for some matches there are no smiles so remove those from the dataframe
    #df.drop(df.index[df['Smiles'] == 'N/A'], inplace=True)
    #df.drop(df.index[df['Smiles'] == ' '], inplace=True)
    # for index, row in df.iterrows():
    #     print(df.at[index,"Smiles"])
    print(df)
    return df