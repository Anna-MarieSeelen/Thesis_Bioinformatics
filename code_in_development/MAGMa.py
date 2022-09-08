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
    connection = sqlite3.connect('db.sqlite')
    cursor = connection.cursor()
    cursor.execute('Create Table if not exists Student (name Text, course Text, roll Integer)')


def main():
    path_to_json_file = argv[1]
    save_json_as_csv(path_to_json_file)

if __name__ == "__main__":
    main()