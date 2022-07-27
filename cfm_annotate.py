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
from sys import argv
import sys
import os
import subprocess
from pathlib import Path

# functions

def annotate_peaks(spectrum_file_name, smiles, identifier, abs_mass_tol=0.01):
    """
    Annotates the MS2 peaks given a smiles and a fragmentation spectrum
    :return:
    """
    # print(df.loc["CCMSLIB00004678842"])
    # cfm-annotate.exe <smiles_or_inchi> <spectrum_file> **<id>** <ppm_mass_tol> <abs_mass_tol> 0.01<param_file> <config_file> <output_file
    file_path_out = Path(r"/lustre/BIF/nobackup/seele006/cfm_annotation_out_{0}".format(identifier))
    out_fn = "cfm_annotation_out_{0}".format(identifier)
    if os.path.exists(file_path_out):
        return file_path_out
    cmd = 'cfm-annotate {0} {1} -abs_mass_tol {2} -output_file {3}' \
            .format(smiles, spectrum_file_name, abs_mass_tol, file_path_out)
    e = subprocess.check_call(cmd, shell=True)
    print("EXIT STATUS AND TYPE", e, type(e))
    print("hi")
    return file_path_out

def main():
    """Main function of this module"""
    # path_to_json_file = argv[1]
    #

    #Make PDF
    # pdf = FPDF()
    # pdf.add_page()
    # pdf.set_font("helvetica", size=10)
    # pdf.image(visualize_mol("N[C@@H](CCCCNC(N)=O)C(O)=O"))
    # pdf.output("output.pdf")
    smiles="OC(C(OC)=C1)=CC=C1C2CC(C3=C(O)C=C(O[C@H]4[C@H](O[C@H]5[C@H](O)[C@H](O)[C@@H](O)[C@H](C)O5)[C@@H](O)[C@H](O)[C@@H](CO)O4)C=C3O2)=O"
    spectrum_file_name='/lustre/BIF/nobackup/seele006/spectrum_files_of_MassQL_matches/spectrum_file_CCMSLIB00006126912.txt'
    identifier="CCMSLIB00006126912"
    annotate_peaks(spectrum_file_name, smiles, identifier)
    return None

if __name__ == "__main__":
    main()