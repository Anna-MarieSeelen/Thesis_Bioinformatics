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

#the older versions of the query maker

def make_MassQL_search2(fragments: list) -> str:
    """
    Takes a list of lists containing fragments and probabilities and returns corresponding MassQL query

    :param fragments: list of lists containing fragments and probabilities belonging to a Mass2Motif
    :return: str, the corresponding MassQL query
    """
    query="QUERY scaninfo(MS2DATA) WHERE POLARITY = Positive " #MS1DATA doesn't give any results, so MS2DATA it is
    fragments=sorted(fragments, key = lambda x: x[1], reverse=True)
    fragments_rel_intensity=list(map(list, fragments)) # to make a copy, but not a reference of fragments
    for i in range(len(fragments_rel_intensity)):
        if i==0:
            fragments_rel_intensity[i][1]=1.0
        else:
            fragments_rel_intensity[i][1]=round(fragments[i][1]/fragments[0][1],3)
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
                        query+="AND MS2NL = {0}:TOLERANCEMZ={1}:INTENSITYMATCH=Y*{2}:INTENSITYMATCHPERCENT={3} ".format(re.search(r'\_(.*)', string).group(1), 0.01,fragment[1],99)
                    else:
                        query+="AND MS2PROD = {0}:TOLERANCEMZ={1}:INTENSITYMATCH=Y*{2}:INTENSITYMATCHPERCENT={3} ".format(re.search(r'\_(.*)', string).group(1), 0.01,fragment[1],99)
                else:
                    pass
    return query

def main() -> None:
    """Main function of this module"""
    fragment_list=argv[1]
    MassQL_query = make_MassQL_search2(fragment_list)
    print(MassQL_query)

if __name__ == "__main__":
    main()