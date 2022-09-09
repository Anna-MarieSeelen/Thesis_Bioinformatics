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
import sqlite3
from sys import argv

#functions
def convert_json_to_sqlite(path_to_json_file):
    connection = sqlite3.connect('/lustre/BIF/nobackup/seele006/MassQL_speclib/HMDB.sqlite')
    cursor = connection.cursor()
    cursor.execute('Create Table if not exists Spectrum (spectrum_id Text, peaks_json Text, Ion_mode Text, Compound_Name Text, Precursor_MZ Float, ExactMass Float, Charge Integer, Smiles Text)')

    HMDB=json.load(open(path_to_json_file))
    columns = ['spectrum_id', 'peaks_json', 'Ion_Mode', 'Compound_Name', 'Precursor_MZ', 'ExactMass', 'Charge', 'Smiles']
    #problem: the peak formatting is not in mgf style...
    for row in HMDB:
        keys = tuple(row[c] for c in columns)
        cursor.execute('insert into Spectrum values(?,?,?,?,?,?,?,?)', keys)
        #print(f'{row["spectrum_id"]} data inserted Succesfully')
    connection.commit()
    connection.close()
    return None

def main():
    path_to_json_file=argv[1]
    convert_json_to_sqlite(path_to_json_file)

if __name__ == "__main__":
    main()