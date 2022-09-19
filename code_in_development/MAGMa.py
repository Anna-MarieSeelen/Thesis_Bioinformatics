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

def main():
    path_to_json_file=argv[1]
    convert_json_to_sqlite(path_to_json_file)

if __name__ == "__main__":
    main()