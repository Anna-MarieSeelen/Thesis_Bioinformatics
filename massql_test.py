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
sys.path.insert(1, '/lustre/BIF/nobackup/seele006/MS2Query_search_libraries')
#TODO: fix that it imports the file, because the massql function sucksss
import ALL_GNPS_210409_positive_processed_annotated_CF_NPC_classes.txt
import ntpath

# functions

def try_massql(query, file):
    results_df=msql_engine.process_query(query,file) # this function can only look in the pwd, but this function should look in head
    return results_df

def main():
    """Main function of this module"""
    # step 1: parse the protein sequence from file 1 and file 2 into dict
    head, tail = ntpath.split(argv[1])
    path_to_file_with_GNPS_mass_library_txt = argv[1]
    query=("QUERY scaninfo(MS2DATA) WHERE MS2NL=163")
    file="ALL_GNPS_210409_positive_processed_annotated_CF_NPC_classes.txt"
    try_massql(query, file)


if __name__ == "__main__":
    main()