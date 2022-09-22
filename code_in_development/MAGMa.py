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
import re
from pathlib import Path
import os
import subprocess
import sys

#functions
def convert_json_to_sqlite(path_to_json_file,dir_database):
    connection = sqlite3.connect('{0}'.format(dir_database))
    cursor = connection.cursor()
    try:
        cursor.execute("""CREATE TABLE molecules (id TEXT PRIMARY KEY,
                                             mim INTEGER NOT NULL,
                                             charge INTEGER NOT NULL,
                                             natoms INTEGER NOT NULL,
                                             molblock TEXT,
                                             inchikey TEXT,
                                             smiles TEXT,
                                             molform TEXT,
                                             name TEXT,
                                             reference TEXT,
                                             logp INT)""")
        connection.commit()
        print("HMDB_MAGMa.db created")
    except:
        print ("HMDB_MAGMa.db already exists (or error creating it)")
        exit()

    HMDB=json.load(open(path_to_json_file))
    columns = ['spectrum_id', 'peaks_json', 'Ion_Mode', 'Compound_Name', 'Precursor_MZ', 'ExactMass', 'Charge', 'Smiles']
    #problem: the peak formatting is not in mgf style...
    for row in HMDB:
        #keys = tuple(row[c] for c in columns)
        #cursor.execute('insert into Spectrum values(?,?,?,?,?,?,?,?)', keys)
        cursor.execute('''INSERT INTO molecules (id, mim, charge, natoms, molblock, inchikey,
                                     smiles,molform,name,reference,logp) VALUES (?,?,?,?,?,?,?,?,?,?,?)''', (
            hmdb_id,
            int(mim * 1e6),
            charge,
            int(natoms),
            unicode(molblock),
            unicode(inchikey),
            unicode(smiles),
            unicode(molform),
            unicode(molname, 'utf-8', 'xmlcharrefreplace'),
            unicode(hmdb_id),
            int(logp * 10)))
        memstore[inchikey] = (hmdb_id, hmdb_id, ionized)
    conn.commit()
        #print(f'{row["spectrum_id"]} data inserted Succesfully')

    cursor.execute(
        'CREATE INDEX idx_cover ON molecules (charge,mim,natoms,reference,molform,inchikey,smiles,name,molblock,logp)')
    connection.commit()
    return None

def initialize_db_to_save_results(path_to_store_results_db,path_to_spectrum_file):
    #generate a good database name with the spectrum id and motif in it.
    spectrum_name=re.search(r'(.*)(spectrum.*)(.txt)', path_to_spectrum_file).group(2)
    file_path_out = Path(r"{0}/MAGMa_db_{1}.sqlite".format(path_to_store_results_db, spectrum_name))
    #if you try to annotate something in an existing database it will go wrong, so if the database exists remove it.
    if os.path.exists(file_path_out):
        os.remove(file_path_out)
    cmd = 'magma init_db {0}'\
            .format(file_path_out)
    try:
        e = subprocess.check_call(cmd, shell=True, stdout = subprocess.DEVNULL, stderr = subprocess.STDOUT)
    except subprocess.CalledProcessError:
        pass
    return file_path_out

def add_spectrum_to_be_annotated_into_db(path_to_results_db_file, path_to_spectrum_file, abs_intensity_thres=1000, mz_precision_ppm=80, mz_precision_abs=0.01, spectrum_file_type='mgf', ionisation=1):
    # a = Absolute intensity threshold for storing peaks in database (default: 1000)
    # i = Ionisation mode (default: 1)
    # p= Maximum relative m/z error (ppm) (default: 5)
    # q= Maximum absolute m/z error (Da) (default: 0.001)
    cmd = 'magma read_ms_data -f {0} -i {1} -a {2} -p {3} -q {4} {5} {6}'.format(spectrum_file_type,ionisation,abs_intensity_thres, mz_precision_ppm, mz_precision_abs, path_to_spectrum_file, path_to_results_db_file)
    e = subprocess.check_call(cmd, shell=True,stdout = subprocess.DEVNULL, stderr = subprocess.STDOUT)
    #print("EXIT STATUS AND TYPE", e, type(e))
    return None

def annotate_spectrum(path_to_structures_database,path_to_results_db_file, max_num_break_bonds=10, structure_db="hmdb", ncpus=1):
    # b = Maximum number of bond breaks to generate substructures (default: 3)
    # s= --structure_database {pubchem,kegg,hmdb} Retrieve molecules from structure database (default: )
    # n = Number of parallel cpus to use for annotation (default: 1)
    cmd = "magma annotate -b {0} -s {1} -o {2} -n {3} {4}".format(max_num_break_bonds,structure_db,path_to_structures_database,ncpus, path_to_results_db_file)
    e = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
    # e=e.decode('utf-8')
    # for i in e.split('\n'):
    #     print(i)
    # if you want to do something with the output
    return None

#TODO: massql will return a sqlite database: you should check the identifiers! Definitely look at the HMDB identifier

#you need the spectrum file with the motif and the identifier
#you need the file from pdf script where you see the fragments for each motif


def main():
    #main function of the script
    #step 0: parse input
    path_to_structures_database=argv[1]
    path_to_spectrum_file=argv[2]
    path_to_store_results_db=argv[3]
    #step 1: initialize database to save results
    path_to_results_db_file=initialize_db_to_save_results(path_to_store_results_db, path_to_spectrum_file)
    print(path_to_results_db_file)
    # step 2: add spectrum to be annotated into the results database
    add_spectrum_to_be_annotated_into_db(path_to_results_db_file, path_to_spectrum_file, abs_intensity_thres=1000,
                                         mz_precision_ppm=80, mz_precision_abs=0.01, spectrum_file_type='mgf',
                                         ionisation=1)
    #step 3: annotate spectrum and store output in results database
    annotate_spectrum(path_to_structures_database, path_to_results_db_file, max_num_break_bonds=10, structure_db="hmdb",
                      ncpus=1)

if __name__ == "__main__":
    main()