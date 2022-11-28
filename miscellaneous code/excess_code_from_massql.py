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
# HEREEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
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

# HEREEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
# ALSO MAKE FUNCTION FROM JSON TO MGF

    with open(mgf_file, 'r') as spectra_file:
        spectra_from_file = list(load_from_mgf(spectra_file))
        for spectrum in spectra_from_file:
            #print(spectrum.get("spectrumid"))
            if spectrum.get("spectrumid") == spectrum_id:
                # in this function I get an error but don't know why
                save_as_mgf(spectrum, file_path)
    return None

def read_json(json_file):
    """
    Puts a json file with MS/MS spectra retrieved in a Pandas dataframe.

    :param json_file: str, the path to a file with MS/MS spectra retrieved from GNPS in json format
    :return: pandas Dataframe with one row for each spectrum and the information belonging to the spectrum in separate
    columns.
    """
    with open(json_file, 'r') as f:
        dict=json.load(f)
    df=json_normalize(dict)
    #print(df.columns)
    df.set_index("spectrum_id",inplace=True, drop=True)
    print(df)
    return df

def make_mgf_txt_file_for_spectrum(motif: str, spectrum_id: str, path_to_store_spectrum_files: str,
                                       dict_with_mgf_spectra: dict) -> str:
    """
    Writes a spectrum containing the motif to a text file in mgf format using the identifier given by MassQL

    :param motif: str, the Mass2Motif for which the MassQL query was made and for which the spectrum was found, which
    contains the motif
    :param spectrum_id: str, the spectrum_id of GNPS spectra files formatted like this: CCMSLIB00000425029
    :param path_to_store_spectrum_files: str, folder where all the mgf-formatted text files with spectra will be stored.
    :param dict_with_mgf_spectra: dict, with {spectrum_id:record} where each record is a string containing the mgf-style
    spectrum of one compound
    :return: path of the spectrum file with the HMDB identifier of the spectrum and the motif it is supposed to contain
    in the file name.
    """
    # look for the corresponding HMDB identifier using the spectrum_id in dict of the mgf file,
    # because MAGMa works with the HMDB identifier.
    spectrum_record_mgf=dict_with_mgf_spectra[spectrum_id]
    # in the current (10-2022) structures database of HMDB which is used for MAGMa a longer identifier is used, so the
    # older identifier from GNPS needs to be adjusted
    HMDB_id_of_spectrum = re.search(r'(.*)(HMDB:)(HMDB\d*)(-\d*)(.*)', spectrum_record_mgf).group(3)
    adj_HMDB_id = HMDB_id_of_spectrum[:4] + '00' + HMDB_id_of_spectrum[4:]
    # create the name for the spectrum file with the HMDB identifier and the motif that the spectrum contains according
    # to MassQL
    path_to_spectrum_file = Path(fr"{path_to_store_spectrum_files}/spectrum_{motif}_{adj_HMDB_id}.txt")
    # write the mgf spectrum saved in the dict only to a text file if the text file doesn't exist already, if the text
    # file does exist it is assumed to be the right text file in the right formatting
    if os.path.exists(path_to_spectrum_file):
         return os.path.abspath(path_to_spectrum_file)
    else:
        spectrum_file = open(path_to_spectrum_file, "w")
        spectrum_file.write(dict_with_mgf_spectra[spectrum_id])
        spectrum_file.close()
        return os.path.abspath(path_to_spectrum_file)

def write_path_to_file(path_to_file_with_motifs: str, path_to_spectrum_file: str) -> None:
    """
    Makes or appends a path of a spectrum file to a text file with paths of the spectrum files which contain Mass2Motifs

    :param path_to_file_with_motifs: str, the path where the file with the selected motifs is stored
    :param path_to_spectrum_file: str, path of the spectrum file with the HMDB identifier of the spectrum and the motif
    it is supposed to contain in the file name.
    :return:
    """
    # the paths to the spectrum files that are associated with all motifs according to MassQL will be stored in this
    # document. This document will be used by the MAGMa.py script.
    path, filename = os.path.split(path_to_file_with_motifs)

    path_to_txt_file_with_paths = Path(fr"{path}/paths_to_spectrum_files_from_MassQL.txt")
    # if the path to the spectrum names file already exists you should append to this file
    if os.path.exists(path_to_txt_file_with_paths):
        # however, you should only append a path that is not yet in the file!
        with open(path_to_txt_file_with_paths, 'r') as file:
            if path_to_spectrum_file in file.read():
                file.close()
                return None
            else:
                output_file = open(path_to_txt_file_with_paths, "a")
                output_file.write(f"{path_to_spectrum_file}")
                output_file.write("\n")
                output_file.close()
    # if the path to the spectrum names file doesn't exist make one and write the path to the spectrum file in it.
    else:
        output_file = open(path_to_txt_file_with_paths, "w")
        output_file.write(f"{path_to_spectrum_file}")
        output_file.write("\n")
        output_file.close()
    return None

def delete_files(path_to_store_spectrum_files):
    # This function should be used in main...
    for file in os.scandir(path_to_store_spectrum_files):
        os.remove(file.path)
    return None

# query = ("QUERY scaninfo(MS2DATA) WHERE POLARITY = Positive AND MS2PROD = 85.0250:TOLERANCEMZ=0.01:INTENSITYMATCH=Y:INTENSITYMATCHREFERENCE AND MS2PROD = 68.0275:TOLERANCEMZ=0.01:INTENSITYMATCH=Y*0.186:INTENSITYMATCHPERCENT=99 AND MS2PROD = 97.0250:TOLERANCEMZ=0.01:INTENSITYMATCH=Y*0.156:INTENSITYMATCHPERCENT=99")
# query = ("QUERY scaninfo(MS2DATA) WHERE POLARITY = Positive AND MS2NL = 176.0350:TOLERANCEMZ=0.01:INTENSITYMATCH=Y:INTENSITYMATCHREFERENCE AND MS2PROD = 126.0550:TOLERANCEMZ=0.01:INTENSITYMATCH=Y*0.089:INTENSITYMATCHPERCENT=99 AND MS2PROD = 127.0375:TOLERANCEMZ=0.01:INTENSITYMATCH=Y*0.082:INTENSITYMATCHPERCENT=99")
# query = ("QUERY scaninfo(MS2DATA) WHERE POLARITY = Positive AND MS2NL = 46.0050:TOLERANCEMZ=0.005") #motif gnps_motif_38.m2m