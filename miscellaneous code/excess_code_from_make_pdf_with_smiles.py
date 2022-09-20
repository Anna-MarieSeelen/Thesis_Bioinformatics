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

