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
import json
import csv
from sys import argv
import re

#functions
def save_json_as_csv(json_file):
    with open(json_file) as json_file:
        jsondata = json.load(json_file)

    data_file = open('/lustre/BIF/nobackup/seele006/MassQL_speclib/HMDB.csv', 'w', newline='')
    csv_writer = csv.writer(data_file)

    count = 0
    for data in jsondata:
        if data["InChIKey_inchi"] != "":
            data["Identifier"] = data["spectrum_id"]
            del data["spectrum_id"]
            data["MonoisotopicMass"] = data["ExactMass"]
            del data["ExactMass"]
            data["MolecularFormula"] = data["Formula_smiles"]
            del data["Formula_smiles"]
            data["SMILES"] = data["Smiles"]
            del data["Smiles"]
            data["InChI"]=data["INCHI"]
            del data["INCHI"]
            data["InChIKey1"]=re.search(r'(.*)(-.*)(-.*)', data["InChIKey_inchi"]).group(1)
            data["InChIKey2"] = re.search(r'(.*-)(.*)(-.*)', data["InChIKey_inchi"]).group(2)
            data["InChIKey3"] = re.search(r'(.*-)(.*-)(.*)', data["InChIKey_inchi"]).group(3)
            data["Name"] = data["Compound_Name"]
            del data["Compound_Name"]
            data["InChIKey"] = data["InChIKey_inchi"]
            del data["InChIKey_inchi"]
            ls=["source_file", "InChIKey_smiles", "task", "scan", "ms_level", "library_membership", "spectrum_status", "peaks_json", "splash", "submit_user", "Ion_Source", "Compound_Source", "Instrument", "PI", "Data_Collector", "Adduct", "Scan", "Precursor_MZ", "Charge", "CAS_Number", "Pubmed_ID", "INCHI_AUX", "Library_Class", "SpectrumID", "Ion_Mode", "create_time", "task_id", "user_id", "Formula_inchi", "url", "annotation_history"]
            #HMDB: HMDB04095 - 2361 5 - Methoxytryptamine: remove everything before the space
            data["Name"]=re.search(r'(.*)(-\d* )(.*)', data["Name"]).group(3)
            for i in ls:
                del data[i]
            if count == 0:
                header = data.keys()
                csv_writer.writerow(header)
                count += 1
            csv_writer.writerow(data.values())

    data_file.close()
    return None

def main():
    path_to_json_file = argv[1]
    save_json_as_csv(path_to_json_file)

if __name__ == "__main__":
    main()