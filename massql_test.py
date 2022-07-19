#!/usr/bin/env python3
"""
Author: Anna-Marie Seelen
Studentnumber:1008970
Description: alignment of protein sequences to a reference protein sequence
and returning stats
Usage: python3 *name_of_script* *name_input_file* [*output_file*] [gap_penalty]
    name_of_script:
    name_input_file:
    output_file: (optional)
    gap_penalty: default is 8 (optional)
"""

# import statements
from massql import msql_engine
from sys import argv
import sys
#TODO: fix that it imports the file, because the massql function sucksss
#import ALL_GNPS_210409_positive_processed_annotated_CF_NPC_classes.txt
import ntpath
from pathlib import Path
import pyarrow.feather as feather
import os
import pyteomics
import re
import rdkit.Chem as Chem
from rdkit.Chem.Draw import MolToImage
from fpdf import FPDF

# functions
def find(name, path):
    for root, dirs, files in os.walk(path):
        if name in files:
            return os.path.join(root, name)

def try_massql(query, file):
    #sys.path.insert(1, '/lustre/BIF/nobackup/seele006/MS2Query_search_libraries/')
    # you can search an mgf file, so maybe you can search a mgf file of a mass spectral library, but then I also want
    # the annotation, and you don't get that here...
    results=msql_engine.process_query(query,file)
    print(results)
    return results

def read_mgf(file):
    reader=pyteomics.mgf.read_header(file)
    print(reader)
    # for i in reader.header.items:
    #     print(i)
    # for dic in reader.header.items:
    #     for item in dic.items:
    #         print(item)
        # for i in reader:
        # print(header(i))
    #dic=pyteomics.mgf.IndexedMGF().read_header()
    #dic=pyteomics.mgf.read_header(reader)
    #print(spectrum)
    #return spectrum

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

def from_feather_to_df(results_feather):
    read_df = feather.read_feather(results_feather)
    print(read_df)
    return read_df

def main():
    """Main function of this module"""
    # step 1: parse the protein sequence from file 1 and file 2 into dict
    path=argv[1]
    #name_file="FractionProfiling_RPOS_ToF10_PreCheck_LTR_01_DDA.mzML"
    #print(find(path, name_file))
    #head, tail = ntpath.split(argv[1])
    #path_to_file_with_GNPS_mass_library_txt = argv[1]
    query=("QUERY scaninfo(MS2DATA) WHERE POLARITY = Positive AND MS2NL = 46.0050:TOLERANCEPPM=5")
    file = "FractionProfiling_RPOS_ToF10_PreCheck_LTR_01_DDA.mzML"
    results_feather=try_massql(query, path)
    #from_feather_to_df(results_feather)
    #read_mgf(path)
    parse_input(path)
    pdf = FPDF()
    pdf.add_page()
    pdf.set_font("helvetica", size=10)
    pdf.image(visualize_mol("N[C@@H](CCCCNC(N)=O)C(O)=O"))
    pdf.output("output.pdf")


if __name__ == "__main__":
    main()